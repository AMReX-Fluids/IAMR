//BL_COPYRIGHT_NOTICE

//
// $Id: Diffusion.cpp,v 1.80 1999-03-25 18:58:40 marc Exp $
//

//
// Comment out this line to use diffusion class outside
// the context of NavierStokes and classes derived from it.
//
#define USE_NAVIERSTOKES 1

#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <Misc.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <BoxDomain.H>
#include <ParmParse.H>
#include <ErrorList.H>
#ifdef USE_NAVIERSTOKES
#include <NavierStokes.H>
#endif
#include <Diffusion.H>
#include <MultiGrid.H>
#include <CGSolver.H>

#include <DIFFUSION_F.H>
#include <VISCOPERATOR_F.H>

#ifdef BL_USE_NEW_HFILES
#include <cfloat>
#else
#include <float.h>
#endif

#if defined(BL_OSF1)
#if defined(BL_USE_DOUBLE)
const Real BL_BOGUS      = DBL_QNAN;
#else
const Real BL_BOGUS      = FLT_QNAN;
#endif
#else
const Real BL_BOGUS      = 1.e30;
#endif

const Real BL_SAFE_BOGUS = -666.e30;

//
// Include files for tensor solve.
//
#include <DivVis.H>
#include <LO_BCTYPES.H>
#include <MCMultiGrid.H>
#include <MCCGSolver.H>
#include <ViscBndryTensor.H>

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define GEOM_GROW 1

Array<int>  Diffusion::is_diffusive;
Array<Real> Diffusion::visc_coef;
Real      Diffusion::visc_tol = 1.0e-10;      // tolerance for viscous solve
Real      Diffusion::visc_abs_tol = 1.0e-10;  // absolute tol. for visc solve

int  Diffusion::first               = 1;
int  Diffusion::do_reflux           = 1;
int  Diffusion::use_cg_solve        = 0;
int  Diffusion::use_tensor_cg_solve = 0;
bool Diffusion::use_mg_precond_flag = false;
int  Diffusion::verbose             = 0;
int  Diffusion::max_order           = 2;
int  Diffusion::tensor_max_order    = 2;
int  Diffusion::scale_abec          = 0;
int  Diffusion::est_visc_mag        = 1;
int  Diffusion::Lphi_in_abs_tol     = 0;

Array<Real> Diffusion::typical_vals;
const Real typical_vals_DEF = 1.0;

Diffusion::Diffusion (Amr*               Parent,
                      AmrLevel*          Caller,
                      Diffusion*         Coarser,
                      int                num_state,
                      FluxRegister*      Viscflux_reg,
                      MultiFab&          Volume,
                      MultiFab*          Area,
                      const Array<int>&  _is_diffusive,
                      const Array<Real>& _visc_coef)
    :
    parent(Parent),
    caller(Caller),
    grids(caller->boxArray()),
    level(caller->Level()),
    coarser(Coarser),
    finer(0),
    NUM_STATE(num_state),
    viscflux_reg(Viscflux_reg),
    volume(Volume),
    area(Area)
{
    if (first)
    {
        first = 0;

        ParmParse ppdiff("diffuse");

        ppdiff.query("v",verbose);

        ppdiff.query("use_cg_solve",use_cg_solve);
        ppdiff.query("use_tensor_cg_solve",use_tensor_cg_solve);
        int use_mg_precond = 0;
        ppdiff.query("use_mg_precond",use_mg_precond);
        use_mg_precond_flag = (use_mg_precond ? true : false);
        ppdiff.query("max_order",max_order);
        ppdiff.query("tensor_max_order",tensor_max_order);
        ppdiff.query("scale_abec",scale_abec);
        ppdiff.query("est_visc_mag",est_visc_mag);
        ppdiff.query("Lphi_in_abs_tol",Lphi_in_abs_tol);

        ParmParse pp("ns");

        pp.query("do_reflux",do_reflux);
        do_reflux = (do_reflux ? 1 : 0);

        pp.query("visc_tol",visc_tol);
        pp.query("visc_abs_tol",visc_abs_tol);

        const int n_visc = _visc_coef.length();
        const int n_diff = _is_diffusive.length();

        if (n_diff < NUM_STATE || n_visc < NUM_STATE)
            BoxLib::Abort("Diffusion::Diffusion(): is_diffusive and/or visc_coef arrays are not long enough");

        visc_coef.resize(NUM_STATE);
        is_diffusive.resize(NUM_STATE);
        for (int i = 0; i < NUM_STATE; i++)
        {
            is_diffusive[i] = _is_diffusive[i];
            visc_coef[i] = _visc_coef[i];
        }
        //
        // Read in typical state sizes.
        //
        typical_vals.resize(NUM_STATE);
        for (int i = 0; i < NUM_STATE; i++)
        {
            typical_vals[i] = typical_vals_DEF;
        }

        int n_typical_vals = Min(NUM_STATE, ppdiff.countval("typical_vals"));
        if (n_typical_vals > 0)
        {
            ppdiff.queryarr("typical_vals",typical_vals,0,n_typical_vals);
        }

        echo_settings();
    }

    if (level > 0)
    {
        crse_ratio = parent->refRatio(level-1);
        coarser->finer = this;
    }
}

Diffusion::~Diffusion () {}

void
Diffusion::echo_settings () const
{
    //
    // Print out my settings.
    //
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "Diffusion settings...\n";
        cout << "  From diffuse:\n";
        cout << "   use_cg_solve =        " << use_cg_solve << '\n';
        cout << "   use_tensor_cg_solve = " << use_tensor_cg_solve << '\n';
        cout << "   use_mg_precond_flag = " << use_mg_precond_flag << '\n';
        cout << "   max_order =           " << max_order << '\n';
        cout << "   tensor_max_order =    " << tensor_max_order << '\n';
        cout << "   scale_abec =          " << scale_abec << '\n';
        cout << "   est_visc_mag =        " << est_visc_mag << '\n';
        cout << "   Lphi_in_abs_tol =     " << Lphi_in_abs_tol << '\n';
    
        cout << "   typical_vals =";
        for (int i = 0; i <NUM_STATE; i++)
            cout << "  " << typical_vals[i];

        cout << "\n\n  From ns:\n";
        cout << "   do_reflux =           " << do_reflux << '\n';
        cout << "   visc_tol =            " << visc_tol << '\n';
        cout << "   visc_abs_tol =        " << visc_abs_tol << '\n';
    
        cout << "   is_diffusive =";
        for (int i =0; i < NUM_STATE; i++)
            cout << "  " << is_diffusive[i];
    
        cout << "\n   visc_coef =";
        for (int i = 0; i < NUM_STATE; i++)
            cout << "  " << visc_coef[i];

        cout << '\n';
    }
}

Real
Diffusion::get_scaled_abs_tol (int                    sigma,
                               const MultiFab*        rhs,
                               Real                   a,
                               Real                   b,
                               const MultiFab*        alpha,
                               const MultiFab* const* betan,
                               const MultiFab* const* betanp1,
                               Real                   reduction) const
{
    //
    // Get norm of rhs.
    //
    Real norm_est, norm_rhs = 0;

    if (rhs != 0)
    {
        assert(grids == rhs->boxArray());

        for (ConstMultiFabIterator Rhsmfi(*rhs); Rhsmfi.isValid(); ++Rhsmfi)
        {
            norm_rhs = Max(norm_rhs,Rhsmfi().norm(0));
        }
        ParallelDescriptor::ReduceRealMax(norm_rhs);
    }

    if (!Lphi_in_abs_tol)
    {
        norm_est = norm_rhs;
    }
    else
    {
        //
        // Approximate (||A||.||x|| + ||rhs||)*reduction
        //
        const Real* dx = caller->Geom().CellSize();
        //
        // Are there spatially varying coefficients?
        //
        int allthere_n, allthere_np1;
        checkBeta(betan,   allthere_n);
        checkBeta(betanp1, allthere_np1);
    
        int do_const_visc = !est_visc_mag || (!allthere_n && !allthere_np1); 
        Real norm_b_div_beta_grad = 0;
        for (int d = 0; d < BL_SPACEDIM; d++)
        {
            Real norm_visc = 0;

            if (do_const_visc)
            {
                norm_visc = Abs( visc_coef[sigma] );
            }
            else
            {
                if (allthere_n)
                {
                    assert(grids == (*betan[d]).boxArray());

                    for (ConstMultiFabIterator Betanmfi(*betan[d]); 
                         Betanmfi.isValid(); ++Betanmfi)
                    {
                        norm_visc = Max(norm_visc,Betanmfi().norm(0));
                    }
                    ParallelDescriptor::ReduceRealMax(norm_visc);

                }
                if (allthere_np1)
                {
                    assert(grids == (*betanp1[d]).boxArray());
                    for (ConstMultiFabIterator Betanp1mfi(*betanp1[d]); 
                         Betanp1mfi.isValid(); ++Betanp1mfi)
                    {
                        norm_visc = Max(norm_visc,Betanp1mfi().norm(0));
                    }
                    ParallelDescriptor::ReduceRealMax(norm_visc);
                }
            }
    
            norm_b_div_beta_grad = Max(norm_b_div_beta_grad,
                                       Abs(b)*norm_visc/(dx[d]*dx[d]));
        }
        //
        // Get norm of alpha * a (if alpha=0, leave this zero).
        //
        Real norm_a_alpha = 0;

        if (alpha != 0)
        {
            assert(grids == alpha->boxArray());

            for (ConstMultiFabIterator Alphamfi(*alpha); Alphamfi.isValid(); ++Alphamfi)
            {
                norm_a_alpha = Max(norm_a_alpha,a * Alphamfi().norm(0));
            }
            ParallelDescriptor::ReduceRealMax(norm_a_alpha);
        }
        //
        // Get norm_A.
        //
        Real norm_A = Max(norm_a_alpha, norm_b_div_beta_grad);
        //
        // Get norm of soln.
        //
        Real norm_x = typical_vals[sigma];

        norm_est = Max( norm_A * norm_x, norm_rhs );
    }
    
    assert(norm_est >= 0);
    
    return norm_est*reduction;
}

void
Diffusion::diffuse_scalar (Real                   dt,
                           int                    sigma,
                           Real                   be_cn_theta,
                           const MultiFab*        rho_half,
                           int                    rho_flag,
                           MultiFab* const*       flux,
                           int                    dataComp,
                           MultiFab*              delta_rhs, 
                           const MultiFab*        alpha, 
                           const MultiFab* const* betan, 
                           const MultiFab* const* betanp1,
                           const SolveMode&       solve_mode)
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
    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "... diffusing scalar: "
             << caller->get_desc_lst()[State_Type].name(sigma) << '\n';

    int allnull, allthere;
    checkBetas(betan, betanp1, allthere, allnull);

    assert(solve_mode==ONEPASS || (delta_rhs && delta_rhs->boxArray()==grids));
    
    const int finest_level = parent->finestLevel();
    //
    // At this point, S_old has bndry at time N, S_new has bndry at time N+1
    //
    MultiFab& S_old = caller->get_old_data(State_Type);
    MultiFab& S_new = caller->get_new_data(State_Type);

    MultiFab Rhs(grids,1,0);
    MultiFab Soln(grids,1,1);
    {
        //
        // Set up Rhs.
        //
        const Real a = 0.0;
        Real b       = -(1.0-be_cn_theta)*dt;
        if (allnull)
            b *= visc_coef[sigma];

        ViscBndry visc_bndry;
        const Real prev_time = caller->get_state_data(State_Type).prevTime();
        ABecLaplacian* visc_op =
            getViscOp(sigma,a,b,prev_time,visc_bndry,rho_half,rho_flag,
                      0,dataComp,betan);
        visc_op->maxOrder(max_order);
        //
        // Copy to single-component multifab, then apply op to rho-scaled state
        //
        MultiFab::Copy(Soln,S_old,sigma,0,1,0);
        if (rho_flag == 2)
            for (MultiFabIterator Smfi(Soln); Smfi.isValid(); ++Smfi)
                Smfi().divide(S_old[Smfi.index()],Smfi.validbox(),Density,0,1);
        visc_op->apply(Rhs,Soln);
        visc_op->compFlux(D_DECL(*flux[0],*flux[1],*flux[2]),Soln);
        for (int i = 0; i < BL_SPACEDIM; ++i)
            (*flux[i]).mult(-b/(dt*caller->Geom().CellSize()[i]));
        delete visc_op;
        //
        // If this is a predictor step, put "explicit" updates passed via S_new
        //  into delta_rhs after scaling by rho_half if reqd, so they dont get lost,
        //  pull it off S_new to avoid double counting
        //   (for rho_flag == 1:
        //       S_new = S_old - dt.(U.Grad(phi)); want Rhs -= rho_half.(U.Grad(phi)),
        //    else
        //       S_new = S_old - dt.Div(U.Phi),   want Rhs -= Div(U.Phi) )
        //
        FArrayBox tmpfab;

        if (solve_mode == PREDICTOR)
        {
            for (MultiFabIterator Smfi(S_new); Smfi.isValid(); ++Smfi)
            {
                const int iGrid = Smfi.index();
                const Box& box  = Smfi.validbox();
                tmpfab.resize(box,1);
                tmpfab.copy(Smfi(),box,sigma,box,0,1);
                tmpfab.minus(S_old[iGrid],box,sigma,0,1);
                Smfi().minus(tmpfab,box,0,sigma,1); // Remove this term from S_new
                tmpfab.mult(1.0/dt,box,0,1);
                if (rho_flag == 1)
                    tmpfab.mult((*rho_half)[iGrid],box,0,0,1);
                if (alpha!=0)
                    tmpfab.mult((*alpha)[iGrid],box,dataComp,0,1);            
                (*delta_rhs)[iGrid].plus(tmpfab,box,0,dataComp,1);
            }
        }
        //
        // Add body sources
        //
        if (delta_rhs != 0)
        {
            for (ConstMultiFabIterator mfi(*delta_rhs); mfi.isValid(); ++mfi)
            {
                const Box& box   = mfi.validbox();
                const int  iGrid = mfi.index();
                tmpfab.resize(box,1);
                tmpfab.copy(mfi(),box,dataComp,box,0,1);
                tmpfab.mult(dt,box,0,1);
                tmpfab.mult(volume[iGrid],box,0,0,1);
                Rhs[iGrid].plus(tmpfab,box,0,0,1);
            }
        }
        //
        // Add hoop stress for x-velocity in r-z coordinates
        // Note: we have to add hoop stress explicitly because the hoop
        // stress which is added through the operator in getViscOp
        // is eliminated by setting a = 0.
        //
#if (BL_SPACEDIM == 2) 
        if (sigma == Xvel && CoordSys::IsRZ())
        {
            for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
            {
                DependentMultiFabIterator volumemfi(Rhsmfi, volume);
                DependentMultiFabIterator S_oldmfi(Rhsmfi, S_old);
                assert(Rhs.box(Rhsmfi.index()) == Rhsmfi.validbox());
                assert(S_old.box(Rhsmfi.index()) == S_oldmfi.validbox());
                assert(volume.box(Rhsmfi.index()) == volumemfi.validbox());
        
                const Box& bx = Rhsmfi.validbox();
                Box sbx       = ::grow(S_oldmfi.validbox(),S_old.nGrow());
                Array<Real> rcen(bx.length(0));
                parent->Geom(level).GetCellLoc(rcen, bx, 0);
                const int*  lo      = bx.loVect();
                const int*  hi      = bx.hiVect();
                const int*  slo     = sbx.loVect();
                const int*  shi     = sbx.hiVect();
                Real*       rhs     = Rhsmfi().dataPtr();
                const Real* sdat    = S_oldmfi().dataPtr(sigma);
                const Real* rcendat = rcen.dataPtr();
                const Real  coeff   = (1.0-be_cn_theta)*visc_coef[sigma]*dt;
                const Real* voli    = volumemfi().dataPtr();
                Box         vbox    = ::grow(volumemfi.validbox(),volume.nGrow());
                const int*  vlo     = vbox.loVect();
                const int*  vhi     = vbox.hiVect();
                FORT_HOOPRHS(rhs, ARLIM(lo), ARLIM(hi), 
                             sdat, ARLIM(slo), ARLIM(shi),
                             rcendat, &coeff, voli, ARLIM(vlo),ARLIM(vhi));
            }
        }
#endif
        //
        // Increment Rhs with S_old*V (or S_old*V*rho_half if rho_flag==1)
        //  (Note: here S_new holds S_old, but also maybe an explicit increment
        //         from advection if solve_mode != PREDICTOR)
        //
        MultiFab::Copy(Soln,S_new,sigma,0,1,0);

        for (MultiFabIterator mfi(Soln); mfi.isValid(); ++mfi)
        {
            const int  iGrid = mfi.index();
            const Box& box   = mfi.validbox();
            mfi().mult(volume[iGrid],box,0,0,1);
            if (rho_flag == 1)
                mfi().mult((*rho_half)[iGrid],box,0,0,1);
            if (alpha!=0)
                mfi().mult((*alpha)[iGrid],box,dataComp,0,1);
            Rhs[iGrid].plus(mfi(),box,0,0,1);
        }
    }
    //
    // Construct viscous operator with bndry data at time N+1.
    //
    const Real a = 1.0;
    Real b       = be_cn_theta*dt;
    if (allnull)
        b *= visc_coef[sigma];

    ViscBndry visc_bndry;
    Real rhsscale = 1.0;
    const Real cur_time = caller->get_state_data(State_Type).curTime();
    ABecLaplacian* visc_op = getViscOp(sigma,a,b,cur_time,visc_bndry,rho_half,
                                       rho_flag,&rhsscale,dataComp,betanp1,alpha);
    Rhs.mult(rhsscale,0,1);
    visc_op->maxOrder(max_order);
    //
    // Make a good guess for Soln
    //
    MultiFab::Copy(Soln,S_new,sigma,0,1,0);
    if (rho_flag == 2)
        for (MultiFabIterator Smfi(Soln); Smfi.isValid(); ++Smfi)
            Smfi().divide(S_new[Smfi.index()],Smfi.validbox(),Density,0,1);
    //
    // Construct solver and call it.
    //
    const Real S_tol     = visc_tol;
    const Real S_tol_abs = get_scaled_abs_tol(sigma, &Rhs, a, b, alpha,
                                              betan, betanp1, visc_abs_tol);
    if (use_cg_solve)
    {
        CGSolver cg(*visc_op,use_mg_precond_flag);
        cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
    else
    {
        MultiGrid mg(*visc_op);
        mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
    //
    // Add in flux contrib from new-time op
    //
    MultiFab** emfSC; // Temporary single-component, edge-based multifab
    allocFluxBoxesLevel(emfSC,0,1);
    visc_op->compFlux(D_DECL(*emfSC[0],*emfSC[1],*emfSC[2]),Soln);
    for (int i = 0; i < BL_SPACEDIM; ++i)
        (*emfSC[i]).mult(b/(dt*caller->Geom().CellSize()[i]));
    delete visc_op;
    for (MultiFabIterator mfi(Soln); mfi.isValid(); ++mfi)
        for (int d = 0; d < BL_SPACEDIM; ++d)
            (*flux[d])[mfi.index()].plus((*emfSC[d])[mfi.index()],0,0,1);
    removeFluxBoxesLevel(emfSC);
    //
    // Copy into state variable at new time, without bc's
    //
    MultiFab::Copy(S_new,Soln,0,sigma,1,0);
    
    if (rho_flag == 2)
        for (MultiFabIterator Smfi(S_new); Smfi.isValid(); ++Smfi)
            Smfi().mult(S_new[Smfi.index()],Smfi.validbox(),Density,sigma,1);
}

void
Diffusion::diffuse_velocity (Real                   dt,
                             Real                   be_cn_theta,
                             const MultiFab*        rho_half,
                             int                    rho_flag,
                             MultiFab*              delta_rhs,
                             const MultiFab* const* betan, 
                             const MultiFab* const* betanp1)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "... diffuse_velocity\n";

    int allnull, allthere;
    checkBetas(betan, betanp1, allthere, allnull);

#ifndef NDEBUG
    for (int d = 0; d < BL_SPACEDIM; ++d)
        assert(allnull ? visc_coef[Xvel+d]>=0 : betan[d]->min(0,0) >= 0.0);
#endif

    if (allnull)
    {
        MultiFab* *fluxSC;
        allocFluxBoxesLevel(fluxSC,0,1);

        for (int sigma = 0; sigma < BL_SPACEDIM; ++sigma)
        {
            const int state_ind = Xvel + sigma;
        
            diffuse_scalar(dt,state_ind,be_cn_theta,rho_half,rho_flag,
                           fluxSC,sigma,delta_rhs);

            if (do_reflux)
            {
                const int finest_level = parent->finestLevel();

                for (int d = 0; d < BL_SPACEDIM; ++d)
                {
                    for (MultiFabIterator fmfi(*fluxSC[d]); fmfi.isValid(); ++fmfi)
                    {
                        if (level < finest_level)
                            finer->viscflux_reg->CrseInit(fmfi(),fmfi().box(),
                                                          d,0,sigma,1,-dt);
                        if (level > 0)
                            viscflux_reg->FineAdd(fmfi(),d,fmfi.index(),0,sigma,1,dt);
                    }
                }
                if (level < finest_level)
                    finer->viscflux_reg->CrseInitFinish();
            }
        }
        removeFluxBoxesLevel(fluxSC);
    }
    else
    {
        diffuse_tensor_velocity(dt,be_cn_theta,rho_half,delta_rhs,betan,betanp1);
    }
}

void
Diffusion::diffuse_tensor_velocity (Real                   dt,
                                    Real                   be_cn_theta,
                                    const MultiFab*        rho_half,
                                    MultiFab*              delta_rhs,
                                    const MultiFab* const* betan, 
                                    const MultiFab* const* betanp1)
{
    const int finest_level = parent->finestLevel();
    //
    // At this point, S_old has bndry at time N S_new contains GRAD(SU).
    //
    MultiFab&  U_old     = caller->get_old_data(State_Type);
    MultiFab&  U_new     = caller->get_new_data(State_Type);
    const Real cur_time  = caller->get_state_data(State_Type).curTime();
    const Real prev_time = caller->get_state_data(State_Type).prevTime();
    const int  dComp     = 0; // FIXME: Start comp for data: dR, betas, should pass in

    int allnull, allthere;
    checkBetas(betan, betanp1, allthere, allnull);
    //
    // U_new now contains the inviscid update of U.
    // This is part of the RHS for the viscous solve.
    //
    MultiFab Rhs(grids,BL_SPACEDIM,0);
    Rhs.setVal(0.0);

    MultiFab** tensorflux_old;
    {
        //
        // Set up Rhs.
        //
        const int soln_old_grow = 1;
        MultiFab Soln_old(grids,BL_SPACEDIM,soln_old_grow);
        const Real a = 0.0;
        Real       b = -(1.0-be_cn_theta)*dt;
        if (allnull)
            b *= visc_coef[Xvel];
        ViscBndryTensor visc_bndry;
        DivVis* tensor_op = getTensorOp(a,b,prev_time,visc_bndry,rho_half,dComp,betan);
        tensor_op->maxOrder(tensor_max_order);
        //
        // Copy to single-component multifab.  Use Soln as a temporary here.
        //
        MultiFab::Copy(Soln_old,U_old,Xvel,0,BL_SPACEDIM,soln_old_grow);
        tensor_op->apply(Rhs,Soln_old);

        if (do_reflux && (level<finest_level || level>0))
        {
            allocFluxBoxesLevel(tensorflux_old,0,BL_SPACEDIM);
            tensor_op->compFlux(*(tensorflux_old[0]),
                                *(tensorflux_old[1]),
#if(BL_SPACEDIM==3)
                                *(tensorflux_old[2]),
#endif
                                Soln_old);
            for (int d = 0; d < BL_SPACEDIM; d++)
                tensorflux_old[d]->mult(-b/(dt*caller->Geom().CellSize()[d]),0);
        }
        delete tensor_op;

        for (int comp = 0; comp < BL_SPACEDIM; comp++)
        {
            int sigma = Xvel + comp;
            //
            // Complete Rhs by adding body sources.
            //
            for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
            {
                DependentMultiFabIterator volumemfi(Rhsmfi, volume);
                DependentMultiFabIterator U_newmfi(Rhsmfi, U_new);
                DependentMultiFabIterator rho_halfmfi(Rhsmfi, (*rho_half));
                assert(grids[Rhsmfi.index()] == Rhsmfi.validbox());
                //
                // Scale inviscid part by volume.
                //
                U_newmfi().mult(volumemfi(),Rhsmfi.validbox(),0,sigma,1);
                //
                // Multiply by density at time nph.
                //
                FArrayBox& Rh = rho_halfmfi();
                U_newmfi().mult(Rh,Rhsmfi.validbox(),0,sigma,1);
                //
                // Add to Rhs which contained operator applied to U_old.
                //
                Rhsmfi().plus(U_newmfi(),Rhsmfi.validbox(),sigma,comp,1);
            }

            if (delta_rhs != 0)
            {
                for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
                {
                    DependentMultiFabIterator volumemfi(Rhsmfi, volume);
                    DependentMultiFabIterator delta_rhsmfi(Rhsmfi,*delta_rhs);
                    assert(grids[Rhsmfi.index()] == Rhsmfi.validbox());
                    delta_rhsmfi().mult(dt,comp+dComp,1);
                    delta_rhsmfi().mult(volumemfi(),Rhsmfi.validbox(),0,comp+dComp,1);
                    Rhsmfi().plus(delta_rhsmfi(),Rhsmfi.validbox(),comp+dComp,comp,1);
                }
            }
        }

#if (BL_SPACEDIM == 2) 
        if (CoordSys::IsRZ())
        {
            int fort_xvel_comp = Xvel+1;

            for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
            {
                DependentMultiFabIterator volumemfi(Rhsmfi, volume);
                DependentMultiFabIterator U_oldmfi(Rhsmfi, U_old);
                DependentMultiFabIterator betanp10mfi(Rhsmfi,*betanp1[0]);
                DependentMultiFabIterator betanp11mfi(Rhsmfi,*betanp1[1]);

                assert(Rhs.box(Rhsmfi.index()) == Rhsmfi.validbox());
                assert(U_old.box(Rhsmfi.index()) == U_oldmfi.validbox());
                assert(volume.box(Rhsmfi.index()) == volumemfi.validbox());

                const Box& bx = Rhsmfi.validbox();
                Box sbx       = ::grow(U_oldmfi.validbox(),U_old.nGrow());
                Array<Real> rcen(bx.length(0));
                parent->Geom(level).GetCellLoc(rcen, bx, 0);
                const int*       lo        = bx.loVect();
                const int*       hi        = bx.hiVect();
                const int*       slo       = sbx.loVect();
                const int*       shi       = sbx.hiVect();
                Real*            rhs       = Rhsmfi().dataPtr();
                const Real*      sdat      = U_oldmfi().dataPtr(Xvel);
                const Real*      rcendat   = rcen.dataPtr();
                const Real       coeff     = (1.0-be_cn_theta)*dt;
                const Real*      voli      = volumemfi().dataPtr();
                Box              vbox      = ::grow(volumemfi.validbox(),volume.nGrow());
                const int*       vlo       = vbox.loVect();
                const int*       vhi       = vbox.hiVect();
                const FArrayBox& betax     = betanp10mfi();
                const int*       betax_lo  = betax.loVect();
                const int*       betax_hi  = betax.hiVect();
                const Real*      betax_dat = betax.dataPtr(dComp);
                const FArrayBox& betay     = betanp11mfi();
                const int*       betay_lo  = betay.loVect();
                const int*       betay_hi  = betay.hiVect();
                const Real*      betay_dat = betay.dataPtr(dComp);

                FORT_TENSOR_HOOPRHS(&fort_xvel_comp, rhs, ARLIM(lo), ARLIM(hi), 
                                    sdat, ARLIM(slo), ARLIM(shi),
                                    rcendat, &coeff, 
                                    voli, ARLIM(vlo), ARLIM(vhi),
                                    betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                                    betay_dat,ARLIM(betay_lo),ARLIM(betay_hi));
            }
        }
#endif
    }
    //
    // I am using a ghost cell in Soln even though Bill does not.
    //
    const int soln_grow = 1;
    MultiFab Soln(grids,BL_SPACEDIM,soln_grow);
    Soln.setVal(0.0);
    //
    // Compute guess of solution.
    //
    if (level == 0)
    {
        MultiFab::Copy(Soln,U_old,Xvel,0,BL_SPACEDIM,soln_grow);
    }
    else
    {
        caller->FillCoarsePatch(Soln,0,cur_time,State_Type,Xvel,BL_SPACEDIM);
    }
    //
    // Copy guess into U_new.
    //
    // The new-time operator is initialized with a "guess" for the new-time
    // state.  We intentionally initialize the grow cells with a bogus
    // value to emphasize that the values are not to be considered "valid"
    // (we shouldn't specify any grow cell information), but rather are to
    // filled by the "physics bc's, etc" in the problem-dependent code.  In
    // the course of this filling (typically while generating/filling the
    // BndryData object for the solvers), StateData::filcc is called to get
    // physical bc's.  Here 'something computable' has to already exist in
    // the grow cells (even though filcc ultimately will fill the corner
    // correctly, if applicable).  This is apparently where the
    // `something computable' is to be set.
    //
    int n_comp  = BL_SPACEDIM;
    int n_ghost = 1;
    U_new.setVal(BL_SAFE_BOGUS,Xvel,n_comp,n_ghost);
    n_ghost = 0;
    U_new.copy(Soln,0,Xvel,n_comp);
    //
    // Construct viscous operator with bndry data at time N+1.
    //
    const Real a = 1.0;
    Real       b = be_cn_theta*dt;
    if (allnull)
        b *= visc_coef[Xvel];
       
    ViscBndryTensor visc_bndry;
    DivVis* tensor_op = getTensorOp(a,b,cur_time,visc_bndry,rho_half,dComp,betanp1);
    tensor_op->maxOrder(tensor_max_order);
    const MultiFab* alpha = &(tensor_op->aCoefficients());
    //
    // Construct solver and call it.
    //
    const Real S_tol     = visc_tol;
    const Real S_tol_abs = get_scaled_abs_tol(Xvel, &Rhs, a, b, alpha, betan,
                                              betanp1, visc_abs_tol);
    if (use_tensor_cg_solve)
    {
        const int use_mg_pre = 0;
        MCCGSolver cg(*tensor_op,use_mg_pre);
        cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
    else
    {
        MCMultiGrid mg(*tensor_op);
        mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }

    int visc_op_lev = 0;
    tensor_op->applyBC(Soln,visc_op_lev); // This may not be needed.
    //
    // Copy into state variable at new time.
    //
    n_ghost = soln_grow;
    MultiFab::Copy(U_new,Soln,0,Xvel,n_comp,n_ghost);
    //
    // Modify diffusive fluxes here.
    //
    if (do_reflux && (level<finest_level || level>0))
    {
        MultiFab** tensorflux;
        allocFluxBoxesLevel(tensorflux,0,BL_SPACEDIM);
        tensor_op->compFlux(*(tensorflux[0]),
                            *(tensorflux[1]),
#if(BL_SPACEDIM==3)
                            *(tensorflux[2]),
#endif
                            Soln);
        for (int d = 0; d < BL_SPACEDIM; d++)
        {
            tensorflux[d]->mult(b/(dt*caller->Geom().CellSize()[d]),0);
            tensorflux[d]->plus(*(tensorflux_old[d]),0,BL_SPACEDIM,0);
        }       
        removeFluxBoxesLevel(tensorflux_old);
    
        FArrayBox xflux, yflux, zflux;

        for (int sigma = Xvel; sigma < BL_SPACEDIM+Xvel; sigma++)
        {

            for (MultiFabIterator tensorflux0mfi(*(tensorflux[0]));
                 tensorflux0mfi.isValid(); ++tensorflux0mfi)
            {
                DependentMultiFabIterator tensorflux1mfi(tensorflux0mfi,
                                                         *(tensorflux[1]));
#if (BL_SPACEDIM == 3)
                DependentMultiFabIterator tensorflux2mfi(tensorflux0mfi,
                                                         *(tensorflux[2]));
#endif

                const int  i   = tensorflux0mfi.index();
                const Box& grd = ::enclosedCells(tensorflux0mfi.validbox());

                assert(grd == grids[tensorflux0mfi.index()]);

                const int* lo  = grd.loVect();
                const int* hi  = grd.hiVect();

                Box xflux_bx(grd);
                xflux_bx.surroundingNodes(0);
                xflux.resize(xflux_bx,1);
                xflux.copy(tensorflux0mfi(),sigma,0,1);

                Box yflux_bx(grd);
                yflux_bx.surroundingNodes(1);
                yflux.resize(yflux_bx,1);
                yflux.copy(tensorflux1mfi(),sigma,0,1);
#if (BL_SPACEDIM == 3)
                Box zflux_bx(grd);
                zflux_bx.surroundingNodes(2);
                zflux.resize(zflux_bx,1);
                zflux.copy(tensorflux2mfi(),sigma,0,1);
#endif
                if (level < finest_level)
                {
                    FluxRegister& fr = *finer->viscflux_reg;
                    fr.CrseInit(xflux,xflux_bx,0,0,sigma,1,-dt);
                    fr.CrseInit(yflux,yflux_bx,1,0,sigma,1,-dt);
#if (BL_SPACEDIM == 3)
                    fr.CrseInit(zflux,zflux_bx,2,0,sigma,1,-dt);
#endif
                }
                if (level > 0)
                {
                    viscflux_reg->FineAdd(xflux,0,i,0,sigma,1,dt);
                    viscflux_reg->FineAdd(yflux,1,i,0,sigma,1,dt);
#if (BL_SPACEDIM == 3)
                    viscflux_reg->FineAdd(zflux,2,i,0,sigma,1,dt);
#endif
                }
            }
            if (level < finest_level)
            {
                finer->viscflux_reg->CrseInitFinish();
            }
        }
        removeFluxBoxesLevel(tensorflux);
    }

    delete tensor_op;
}

void
Diffusion::diffuse_Vsync (MultiFab*              Vsync,
                          Real                   dt,
                          Real                   be_cn_theta,
                          const MultiFab*        rho_half,
                          int                    rho_flag,
                          const MultiFab* const* beta)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "Diffusion::diffuse_Vsync\n";

    int allnull, allthere;
    checkBeta(beta, allthere, allnull);

#ifndef NDEBUG
    for (int d = 0; d < BL_SPACEDIM; ++d)
        assert(allnull ? visc_coef[Xvel+d]>=0 : beta[d]->min(0,0) >= 0.0);
#endif

    if (allnull)
        diffuse_Vsync_constant_mu(Vsync,dt,be_cn_theta,rho_half,rho_flag);
    else
        diffuse_tensor_Vsync(Vsync,dt,be_cn_theta,rho_half,rho_flag,beta);
    //
    // applyBC has put "incorrect" values in the ghost cells
    // outside external Dirichlet boundaries. Reset these to zero
    // so that syncproject and conservative interpolation works correctly.
    //
    Box domain = ::grow(caller->Geom().Domain(),1);

    for (int n = Xvel; n < Xvel+BL_SPACEDIM; n++)
    {
        const BCRec& velbc = caller->get_desc_lst()[State_Type].getBC(n);

        for (int k = 0; k < BL_SPACEDIM; k++)
        {
            if (velbc.hi(k) == EXT_DIR)
            {
                IntVect smallend = domain.smallEnd();
                smallend.setVal(k,domain.bigEnd(k));
                Box top_strip(smallend,domain.bigEnd(),IntVect::TheCellVector());
                Vsync->setVal(0.0,top_strip,n-Xvel,1,1);
            }
            if (velbc.lo(k) == EXT_DIR)
            {
                IntVect bigend = domain.bigEnd();
                bigend.setVal(k,domain.smallEnd(k));
                Box bottom_strip(domain.smallEnd(),bigend,IntVect::TheCellVector());
                Vsync->setVal(0.0,bottom_strip,n-Xvel,1,1);
            }
        }
    }
}

void
Diffusion::diffuse_Vsync_constant_mu (MultiFab*       Vsync,
                                      Real            dt,
                                      Real            be_cn_theta,
                                      const MultiFab* rho_half,
                                      int             rho_flag)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "Diffusion::diffuse_Vsync\n";

    const Real* dx     = caller->Geom().CellSize();
    const int   IOProc = ParallelDescriptor::IOProcessorNumber();
    //
    // At this point in time we can only do decoupled scalar
    // so we loop over components.
    //
    for (int comp = 0; comp < BL_SPACEDIM; comp++)
    {
        MultiFab Soln(grids,1,1);
        MultiFab Rhs(grids,1,0);

        Soln.setVal(0);
        Rhs.copy(*Vsync,comp,0,1);

        if (verbose)
        {
            Real r_norm = 0.0;
            for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
            {
                r_norm = Max(r_norm,Rhsmfi().norm(0));
            }
            ParallelDescriptor::ReduceRealMax(r_norm,IOProc);

            if (ParallelDescriptor::IOProcessor())
                cout << "Original max of Vsync " << r_norm << '\n';
        }
        //
        // Multiply RHS by volume and density.
        //
        for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
        {
            DependentMultiFabIterator volumemfi(Rhsmfi, volume);
            DependentMultiFabIterator rho_halfmfi(Rhsmfi, (*rho_half));
            Rhsmfi().mult(volumemfi()); 
            Rhsmfi().mult(rho_halfmfi()); 
        }
        //
        // SET UP COEFFICIENTS FOR VISCOUS SOLVER.
        //
        const Real a  = 1.0;
        const Real b  = be_cn_theta*dt*visc_coef[comp];
        Real rhsscale = 1.0;
        ABecLaplacian* visc_op = getViscOp(comp,a,b,rho_half,rho_flag,&rhsscale);
        visc_op->maxOrder(max_order);
        Rhs.mult(rhsscale,0,1);
        //
        // Construct solver and call it.
        //
        const Real S_tol      = visc_tol;
        const MultiFab* alpha = &(visc_op->aCoefficients());

        MultiFab const* betan[BL_SPACEDIM];
        MultiFab const* betanp1[BL_SPACEDIM];

        for (int d = 0; d < BL_SPACEDIM; d++)
        {
            betan[d] = &visc_op->bCoefficients(d);
            betanp1[d] = betan[d];
        }
        const Real S_tol_abs = get_scaled_abs_tol(comp, &Rhs, a, b, alpha,
                                                  betan, betanp1, visc_abs_tol);
        if (use_cg_solve)
        {
            CGSolver cg(*visc_op,use_mg_precond_flag);
            cg.solve(Soln,Rhs,S_tol,S_tol_abs);
        }
        else
        {
            MultiGrid mg(*visc_op);
            mg.solve(Soln,Rhs,S_tol,S_tol_abs);
        }

        int visc_op_lev = 0;
        visc_op->applyBC(Soln,visc_op_lev);

        MultiFab::Copy(*Vsync,Soln,0,comp,1,1);

        if (verbose)
        {
            Real s_norm = 0.0;
            for (MultiFabIterator Solnmfi(Soln); Solnmfi.isValid(); ++Solnmfi)
            {
                s_norm = Max(s_norm,Solnmfi().norm(0));
            }
            ParallelDescriptor::ReduceRealMax(s_norm,IOProc);

            if (ParallelDescriptor::IOProcessor())
                cout << "Final max of Vsync " << s_norm << '\n';
        }

        delete visc_op;

        FArrayBox xflux, yflux;

        if (level > 0)
        {
            for (MultiFabIterator Vsyncmfi(*Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
            {
                DependentMultiFabIterator area0mfi(Vsyncmfi, area[0]);
                DependentMultiFabIterator area1mfi(Vsyncmfi, area[1]);
#if (BL_SPACEDIM == 3)
                DependentMultiFabIterator area2mfi(Vsyncmfi, area[2]);
#endif
                assert(grids[Vsyncmfi.index()] == Vsyncmfi.validbox());

                const int  i      = Vsyncmfi.index();
                const Box& grd    = Vsyncmfi.validbox();
                const int* lo     = grd.loVect();
                const int* hi     = grd.hiVect();
                FArrayBox& u_sync = Vsyncmfi();
                const int* ulo    = u_sync.loVect();
                const int* uhi    = u_sync.hiVect();

                Box xflux_bx(grd);
                xflux_bx.surroundingNodes(0);
                xflux.resize(xflux_bx,1);
                DEF_LIMITS(xflux,xflux_dat,xflux_lo,xflux_hi);

                Box yflux_bx(grd);
                yflux_bx.surroundingNodes(1);
                yflux.resize(yflux_bx,1);
                DEF_LIMITS(yflux,yflux_dat,yflux_lo,yflux_hi);

                FArrayBox& xarea = area0mfi();
                FArrayBox& yarea = area1mfi();

                DEF_CLIMITS(xarea,xarea_dat,xarea_lo,xarea_hi);
                DEF_CLIMITS(yarea,yarea_dat,yarea_lo,yarea_hi);
                //
                // The extra factor of dt comes from the fact that Vsync
                // looks like dV/dt, not just an increment to V.
                //
                Real mult = -be_cn_theta*dt*dt*visc_coef[comp];

#if (BL_SPACEDIM == 2)
                FORT_VISCSYNCFLUX (u_sync.dataPtr(comp), ARLIM(ulo), ARLIM(uhi),
                                   lo,hi,
                                   xflux_dat,ARLIM(xflux_lo),ARLIM(xflux_hi),
                                   yflux_dat,ARLIM(yflux_lo),ARLIM(yflux_hi),
                                   xarea_dat,ARLIM(xarea_lo),ARLIM(xarea_hi),
                                   yarea_dat,ARLIM(yarea_lo),ARLIM(yarea_hi),
                                   dx,&mult);
#endif
#if (BL_SPACEDIM == 3)

                FArrayBox zflux;
                Box zflux_bx(grd);
                zflux_bx.surroundingNodes(2);
                zflux.resize(zflux_bx,1);
                DEF_LIMITS(zflux,zflux_dat,zflux_lo,zflux_hi);

                FArrayBox& zarea = area2mfi();
                DEF_CLIMITS(zarea,zarea_dat,zarea_lo,zarea_hi);

                FORT_VISCSYNCFLUX (u_sync.dataPtr(comp), ARLIM(ulo), ARLIM(uhi),
                                   lo,hi,
                                   xflux_dat,ARLIM(xflux_lo),ARLIM(xflux_hi),
                                   yflux_dat,ARLIM(yflux_lo),ARLIM(yflux_hi),
                                   zflux_dat,ARLIM(zflux_lo),ARLIM(zflux_hi),
                                   xarea_dat,ARLIM(xarea_lo),ARLIM(xarea_hi),
                                   yarea_dat,ARLIM(yarea_lo),ARLIM(yarea_hi),
                                   zarea_dat,ARLIM(zarea_lo),ARLIM(zarea_hi),
                                   dx,&mult);
#endif

                Real one = 1.0;
                viscflux_reg->FineAdd(xflux,0,i,0,comp,1,one);
                viscflux_reg->FineAdd(yflux,1,i,0,comp,1,one);
#if (BL_SPACEDIM == 3)
                viscflux_reg->FineAdd(zflux,2,i,0,comp,1,one);
#endif
            }
        }
    }
}

void
Diffusion::diffuse_tensor_Vsync (MultiFab*              Vsync,
                                 Real                   dt,
                                 Real                   be_cn_theta,
                                 const MultiFab*        rho_half,
                                 int                    rho_flag,
                                 const MultiFab* const* beta)
{
    const int   finest_level = parent->finestLevel();
    const Real* dx           = caller->Geom().CellSize();
    const int   IOProc       = ParallelDescriptor::IOProcessorNumber();

    MultiFab Soln(grids,BL_SPACEDIM,1);
    MultiFab Rhs(grids,BL_SPACEDIM,0);

    Soln.setVal(0);
    Rhs.copy(*Vsync,0,0,BL_SPACEDIM);

    if (verbose)
    {
        Real r_norm = 0.0;
        for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
        {
            r_norm = Max(r_norm,Rhsmfi().norm(0));
        }
        ParallelDescriptor::ReduceRealMax(r_norm,IOProc);

        if (ParallelDescriptor::IOProcessor())
            cout << "Original max of Vsync " << r_norm << '\n';
    }
    //
    // Multiply RHS by volume and density.
    //
    for (int comp = 0; comp < BL_SPACEDIM; comp++)
    {
        for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
        {
            DependentMultiFabIterator volumemfi(Rhsmfi, volume);
            DependentMultiFabIterator rho_halfmfi(Rhsmfi, (*rho_half));
            Rhsmfi().mult(volumemfi(),0,comp,1); 
            Rhsmfi().mult(rho_halfmfi(),0,comp,1); 
        }
    }
    //
    // SET UP COEFFICIENTS FOR VISCOUS SOLVER.
    //
    const Real a = 1.0;
    const Real b = be_cn_theta*dt;

    const int dComp = 0; // FIXME: start comp for betas, should be passed in
    DivVis* tensor_op = getTensorOp(a,b,rho_half,dComp,beta);
    tensor_op->maxOrder(tensor_max_order);
    //
    // Construct solver and call it.
    //
    const Real S_tol      = visc_tol;
    const MultiFab* alpha = &(tensor_op->aCoefficients());
    MultiFab const* betan[BL_SPACEDIM];
    MultiFab const* betanp1[BL_SPACEDIM];
    for (int d = 0; d < BL_SPACEDIM; d++)
    {
        betan[d]   = &tensor_op->bCoefficients(d);
        betanp1[d] = betan[d];
    }
    const Real S_tol_abs = get_scaled_abs_tol(Xvel, &Rhs, a, b, alpha, betan,
                                              betanp1, visc_abs_tol);

    if (use_tensor_cg_solve)
    {
        MCCGSolver cg(*tensor_op,use_mg_precond_flag);
        cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
    else
    {
        MCMultiGrid mg(*tensor_op);
        mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }

    int visc_op_lev = 0;
    tensor_op->applyBC(Soln,visc_op_lev); 

    MultiFab::Copy(*Vsync,Soln,0,0,BL_SPACEDIM,1);

    if (verbose)
    {
        Real s_norm = 0.0;
        for (MultiFabIterator Solnmfi(Soln); Solnmfi.isValid(); ++Solnmfi)
        {
            s_norm = Max(s_norm,Solnmfi().norm(0));
        }
        ParallelDescriptor::ReduceRealMax(s_norm,IOProc);

        if (ParallelDescriptor::IOProcessor())
            cout << "Final max of Vsync " << s_norm << '\n';
    }

    FArrayBox xflux, yflux, zflux;

    if (level > 0)
    {
        MultiFab** tensorflux;
        allocFluxBoxesLevel(tensorflux,0,BL_SPACEDIM);
        tensor_op->compFlux(*(tensorflux[0]),*(tensorflux[1]),
#if(BL_SPACEDIM==3)
                            *(tensorflux[2]),
#endif
                            Soln);
        //
        // The extra factor of dt comes from the fact that Vsync looks
        // like dV/dt, not just an increment to V.
        //
        // This is to remove the dx scaling in the coeffs
        //
        for (int d =0; d <BL_SPACEDIM; d++)
            tensorflux[d]->mult(b/(dt*caller->Geom().CellSize()[d]),0);

        for (int sigma = Xvel; sigma < BL_SPACEDIM+Xvel; sigma++)
        {
            for (MultiFabIterator tensorflux0mfi(*(tensorflux[0]));
                 tensorflux0mfi.isValid(); ++tensorflux0mfi)
            {
                DependentMultiFabIterator tensorflux1mfi(tensorflux0mfi,
                                                         (*(tensorflux[1])));
#if (BL_SPACEDIM == 3)
                DependentMultiFabIterator tensorflux2mfi(tensorflux0mfi,
                                                         (*(tensorflux[2])));
#endif
                const int i    = tensorflux0mfi.index();
                const Box& grd = ::enclosedCells(tensorflux0mfi.validbox());

                assert(grd==grids[tensorflux0mfi.index()]);

                const int* lo  = grd.loVect();
                const int* hi  = grd.hiVect();

                Box xflux_bx(grd);
                xflux_bx.surroundingNodes(0);
                xflux.resize(xflux_bx,1);
                xflux.copy(tensorflux0mfi(),sigma,0,1);

                Box yflux_bx(grd);
                yflux_bx.surroundingNodes(1);
                yflux.resize(yflux_bx,1);
                yflux.copy(tensorflux1mfi(),sigma,0,1); 
#if (BL_SPACEDIM == 3)
                Box zflux_bx(grd);
                zflux_bx.surroundingNodes(2);
                zflux.resize(zflux_bx,1);
                zflux.copy(tensorflux2mfi(),sigma,0,1);
#endif
                viscflux_reg->FineAdd(xflux,0,i,0,sigma,1,dt);
                viscflux_reg->FineAdd(yflux,1,i,0,sigma,1,dt);
#if (BL_SPACEDIM == 3)
                viscflux_reg->FineAdd(zflux,2,i,0,sigma,1,dt);
#endif
            }
        }
        removeFluxBoxesLevel(tensorflux);
    }
    delete tensor_op;
}

void
Diffusion::diffuse_Ssync (MultiFab*              Ssync,
                          int                    sigma,
                          Real                   dt,
                          Real                   be_cn_theta,
                          const MultiFab*        rho_half,
                          int                    rho_flag,
                          MultiFab* const*       flux,
                          int                    dataComp,
                          const MultiFab* const* beta,
                          const MultiFab*        alpha)
{
    const int state_ind = sigma + BL_SPACEDIM;
    const int IOProc    = ParallelDescriptor::IOProcessorNumber();

    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "Diffusion::diffuse_Ssync: "
             << caller->get_desc_lst()[State_Type].name(state_ind) << '\n';

    int allnull, allthere;
    checkBeta(beta, allthere, allnull);

    MultiFab Soln(grids,1,1);
    MultiFab Rhs(grids,1,0);

    Soln.setVal(0);
    MultiFab::Copy(Rhs,*Ssync,sigma,0,1,0);

    if (verbose)
    {
        MultiFab junk(grids,1,0);

        MultiFab::Copy(junk,Rhs,0,0,1,0);

        if (rho_flag == 2)
        {
            MultiFab& S_new = caller->get_new_data(State_Type);
            for (MultiFabIterator jmfi(junk); jmfi.isValid(); ++jmfi)
            {
                DependentMultiFabIterator S_newmfi(jmfi,S_new);
                jmfi().divide(S_newmfi(),jmfi.validbox(),Density,0,1);
            }
        }
        Real r_norm = 0.0;
        for (MultiFabIterator jmfi(junk); jmfi.isValid(); ++jmfi)
        {
            r_norm = Max(r_norm,jmfi().norm(0));
        }
        ParallelDescriptor::ReduceRealMax(r_norm,IOProc);

        if (ParallelDescriptor::IOProcessor())
            cout << "Original max of Ssync " << r_norm << '\n';
    }
    //
    // Compute RHS.
    //
    for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
        DependentMultiFabIterator volumemfi(Rhsmfi, volume);
        Rhsmfi().mult(volumemfi()); 
        if (rho_flag == 1)
        {
            DependentMultiFabIterator rho_halfmfi(Rhsmfi, *rho_half);
            Rhsmfi().mult(rho_halfmfi());
        }
    }
    //
    // SET UP COEFFICIENTS FOR VISCOUS SOLVER.
    //
    const Real a = 1.0;
    Real       b = be_cn_theta*dt;
    if (allnull)
        b *= visc_coef[state_ind];
    Real rhsscale = 1.0;
    ABecLaplacian* visc_op =
        getViscOp(state_ind,a,b,rho_half,rho_flag,&rhsscale,dataComp,beta,alpha);
    visc_op->maxOrder(max_order);
    Rhs.mult(rhsscale,0,1);
    //
    // Construct solver and call it.
    //
    const Real S_tol = visc_tol;
    MultiFab const * betan[BL_SPACEDIM];
    MultiFab const * betanp1[BL_SPACEDIM];
    for (int d = 0; d < BL_SPACEDIM; d++)
    {
        betan[d]   = &visc_op->bCoefficients(d);
        betanp1[d] = betan[d];
    }
    const Real S_tol_abs = get_scaled_abs_tol(state_ind, &Rhs, a, b,
                                              alpha, betan, betanp1, visc_abs_tol);
    if (use_cg_solve)
    {
        CGSolver cg(*visc_op,use_mg_precond_flag);
        cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }
    else
    {
        MultiGrid mg(*visc_op);
        mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }

    int flux_allthere, flux_allnull;
    checkBeta(flux, flux_allthere, flux_allnull);
    if (flux_allthere)
    {
        visc_op->compFlux(D_DECL(*flux[0],*flux[1],*flux[2]),Soln);
        for (int i = 0; i < BL_SPACEDIM; ++i)
            (*flux[i]).mult(b/(dt*caller->Geom().CellSize()[i]),0);
    }

    MultiFab::Copy(*Ssync,Soln,0,sigma,1,0);
    
    if (verbose)
    {
        Real s_norm = 0.0;
        for (MultiFabIterator Solnmfi(Soln); Solnmfi.isValid(); ++Solnmfi)
        {
            s_norm = Max(s_norm,Solnmfi().norm(0));
        }
        ParallelDescriptor::ReduceRealMax(s_norm,IOProc);

        if (ParallelDescriptor::IOProcessor())
            cout << "Final max of Ssync " << s_norm << '\n';
    }
    
    if (rho_flag == 2)
    {
        MultiFab& S_new = caller->get_new_data(State_Type);

        for (MultiFabIterator Ssyncmfi(*Ssync); Ssyncmfi.isValid(); ++Ssyncmfi)
        {
            DependentMultiFabIterator S_newmfi(Ssyncmfi,S_new);
            Ssyncmfi().mult(S_newmfi(),Ssyncmfi.validbox(),Density,sigma,1);
        }
    }
    
    delete visc_op;
}

DivVis*
Diffusion::getTensorOp (Real                   a,
                        Real                   b,
                        Real                   time,
                        ViscBndryTensor&       visc_bndry,
                        const MultiFab*        rho_half,
                        int                    dataComp,
                        const MultiFab* const* beta)
{
    int allthere;
    checkBeta(beta, allthere);

    const Real* dx = caller->Geom().CellSize();

    getTensorBndryData(visc_bndry,time);

    DivVis* tensor_op = new DivVis(visc_bndry,dx);
    tensor_op->maxOrder(tensor_max_order);

    int isrz   = CoordSys::IsRZ();
    const int nghost = 1; // Just like Bill.
    //
    // alpha should be the same size as volume.
    //
    const int nCompAlpha = BL_SPACEDIM == 2  ?  2  :  1;
    MultiFab alpha(grids,nCompAlpha,nghost);
    alpha.setVal(0.0,nghost);

    if (a != 0.0)
    {
        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator volumemfi(alphamfi, volume);
            DependentMultiFabIterator rho_halfmfi(alphamfi, (*rho_half));
            DependentMultiFabIterator beta0mfi(alphamfi, (*beta[0]));
            DependentMultiFabIterator beta1mfi(alphamfi, (*beta[1]));
#if (BL_SPACEDIM==3)
            DependentMultiFabIterator beta2mfi(alphamfi, (*beta[2]));
#endif
            assert(alpha.box(alphamfi.index()) == alphamfi.validbox());
            assert(volume.box(alphamfi.index()) == volumemfi.validbox());

            const Box& bx       = alphamfi.validbox();
            Array<Real> rcen(bx.length(0));
            parent->Geom(level).GetCellLoc(rcen, bx, 0);
            const int*  lo        = bx.loVect();
            const int*  hi        = bx.hiVect();
            Real*       alpha_dat = alphamfi().dataPtr();
            Box         abx       = ::grow(bx,alpha.nGrow());
            const int*  alo       = abx.loVect();
            const int*  ahi       = abx.hiVect();
            const Real* rcendat   = rcen.dataPtr();
            const Real* voli      = volumemfi().dataPtr();
            Box         vbox      = ::grow(volumemfi.validbox(),volume.nGrow());
            const int*  vlo       = vbox.loVect();
            const int*  vhi       = vbox.hiVect();

            FArrayBox& Rh = rho_halfmfi();
            DEF_LIMITS(Rh,rho_dat,rlo,rhi);

            FArrayBox&  betax     = beta0mfi();
            const Real* betax_dat = betax.dataPtr(dataComp);
            const int*  betax_lo  = betax.loVect();
            const int*  betax_hi  = betax.hiVect();

            FArrayBox&  betay     = beta1mfi();
            const Real* betay_dat = betay.dataPtr(dataComp);
            const int*  betay_lo  = betay.loVect();
            const int*  betay_hi  = betay.hiVect();

#if (BL_SPACEDIM == 3)
            FArrayBox&  betaz     = beta2mfi();
            const Real* betaz_dat = betaz.dataPtr(dataComp);
            const int*  betaz_lo  = betaz.loVect();
            const int*  betaz_hi  = betaz.hiVect();
#endif

            FORT_SET_TENSOR_ALPHA(alpha_dat, ARLIM(alo), ARLIM(ahi),
                                  lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                                  voli, ARLIM(vlo), ARLIM(vhi),
                                  rho_dat,ARLIM(rlo),ARLIM(rhi),
                                  betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                                  betay_dat,ARLIM(betay_lo),ARLIM(betay_hi),
#if (BL_SPACEDIM == 3)
                                  betaz_dat,ARLIM(betaz_lo),ARLIM(betaz_hi),
#endif
                                  &isrz);
        }
    }
    tensor_op->setScalars(a,b);
    tensor_op->aCoefficients(alpha);

    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        MultiFab bcoeffs(area[n].boxArray(),1,nghost);
        bcoeffs.setBndry(0);
        MultiFab::Copy(bcoeffs,area[n],0,0,1,nghost);
        for (MultiFabIterator bcoeffsmfi(bcoeffs); bcoeffsmfi.isValid(); ++bcoeffsmfi)
        {
            bcoeffsmfi().mult(dx[n]);
            bcoeffsmfi().mult((*beta[n])[bcoeffsmfi.index()],dataComp,0,1);
        }
        tensor_op->bCoefficients(bcoeffs,n);
    }

    return tensor_op;
}

DivVis*
Diffusion::getTensorOp (Real                   a,
                        Real                   b,
                        const MultiFab*        rho_half,
                        int                    dataComp,
                        const MultiFab* const* beta)
{
    int allthere = beta != 0;
    if (allthere)
    {
        for (int d = 0; d < BL_SPACEDIM; d++)
        {
            allthere = allthere && beta[d] != 0;
        }
    }
    if (!allthere)
        BoxLib::Abort("Diffusion::getTensorOp(): all betas must allocated all 0 or all non-0");

    const Real* dx     = caller->Geom().CellSize();
    const Box&  domain = caller->Geom().Domain();

    const int nDer = MCLinOp::bcComponentsNeeded();

    Array<BCRec> bcarray(nDer,BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
                                    D_DECL(EXT_DIR,EXT_DIR,EXT_DIR)));

    for (int id = 0; id < BL_SPACEDIM; id++)
    {
        bcarray[id] = caller->get_desc_lst()[State_Type].getBC(Xvel+id);
    }

    IntVect ref_ratio = level > 0 ? parent->refRatio(level-1) : IntVect::TheUnitVector();

    ViscBndryTensor bndry;
    //bndry.define(grids,2*BL_SPACEDIM,caller->Geom());
    bndry.define(grids,nDer,caller->Geom());
    bndry.setHomogValues(bcarray, ref_ratio[0]);
    DivVis* tensor_op = new DivVis(bndry,dx);
    tensor_op->maxOrder(tensor_max_order);

    int isrz   = CoordSys::IsRZ();
    const int nghost = 1; // Just like Bill.
    //
    // alpha should be the same size as volume.
    //
    const int nCompAlpha = BL_SPACEDIM == 2  ?  2  :  1;
    MultiFab alpha(grids,nCompAlpha,nghost);
    alpha.setVal(0.0);

    if (a != 0.0)
    {
        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator volumemfi(alphamfi, volume);
            DependentMultiFabIterator rho_halfmfi(alphamfi, (*rho_half));
            DependentMultiFabIterator beta0mfi(alphamfi, (*beta[0]));
            DependentMultiFabIterator beta1mfi(alphamfi, (*beta[1]));
#if (BL_SPACEDIM==3)
            DependentMultiFabIterator beta2mfi(alphamfi, (*beta[2]));
#endif

            assert(alpha.box(alphamfi.index()) == alphamfi.validbox());
            assert(volume.box(alphamfi.index()) == volumemfi.validbox());

            const Box& bx         = alphamfi.validbox();
            Array<Real> rcen(bx.length(0));
            parent->Geom(level).GetCellLoc(rcen, bx, 0);
            const int*  lo        = bx.loVect();
            const int*  hi        = bx.hiVect();
            Real*       alpha_dat = alphamfi().dataPtr();
            Box         abx       = ::grow(bx,alpha.nGrow());
            const int*  alo       = abx.loVect();
            const int*  ahi       = abx.hiVect();
            const Real* rcendat   = rcen.dataPtr();
            const Real* voli      = volumemfi().dataPtr();
            Box         vbox      = ::grow(volumemfi.validbox(),volume.nGrow());
            const int*  vlo       = vbox.loVect();
            const int*  vhi       = vbox.hiVect();

            FArrayBox& Rh = rho_halfmfi();
            DEF_LIMITS(Rh,rho_dat,rlo,rhi);

            FArrayBox&  betax     = beta0mfi();
            const Real* betax_dat = betax.dataPtr(dataComp);
            const int*  betax_lo  = betax.loVect();
            const int*  betax_hi  = betax.hiVect();
            FArrayBox&  betay     = beta1mfi();
            const Real* betay_dat = betay.dataPtr(dataComp);
            const int*  betay_lo  = betay.loVect();
            const int*  betay_hi  = betay.hiVect();

#if (BL_SPACEDIM == 3)
            FArrayBox&  betaz     = beta2mfi();
            const Real* betaz_dat = betaz.dataPtr(dataComp);
            const int*  betaz_lo  = betaz.loVect();
            const int*  betaz_hi  = betaz.hiVect();
#endif

            FORT_SET_TENSOR_ALPHA(alpha_dat, ARLIM(alo), ARLIM(ahi),
                                  lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                                  voli, ARLIM(vlo), ARLIM(vhi),
                                  rho_dat,ARLIM(rlo),ARLIM(rhi),
                                  betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                                  betay_dat,ARLIM(betay_lo),ARLIM(betay_hi),
#if (BL_SPACEDIM == 3)
                                  betaz_dat,ARLIM(betaz_lo),ARLIM(betaz_hi),
#endif
                                  &isrz);
        }
    }
    tensor_op->setScalars(a,b);
    tensor_op->aCoefficients(alpha);

    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        MultiFab bcoeffs(area[n].boxArray(),1,nghost);
        bcoeffs.setBndry(0);
        MultiFab::Copy(bcoeffs,area[n],0,0,1,nghost);
        for (MultiFabIterator bcoeffsmfi(bcoeffs); bcoeffsmfi.isValid(); ++bcoeffsmfi)
        {
            bcoeffsmfi().mult(dx[n]);
            bcoeffsmfi().mult((*beta[n])[bcoeffsmfi.index()],dataComp,0,1);
        }
        tensor_op->bCoefficients(bcoeffs,n);
    }

    return tensor_op;
}

ABecLaplacian*
Diffusion::getViscOp (int                    comp,
                      Real                   a,
                      Real                   b,
                      Real                   time,
                      ViscBndry&             visc_bndry,
                      const MultiFab*        rho_half,
                      int                    rho_flag, 
                      Real*                  rhsscale,
                      int                    dataComp,
                      const MultiFab* const* beta,
                      const MultiFab*        alpha_in)
{
    int allnull, allthere;
    checkBeta(beta, allthere, allnull);

    const Real* dx = caller->Geom().CellSize();

    getBndryData(visc_bndry,comp,1,time,rho_flag);

    ABecLaplacian* visc_op = new ABecLaplacian(visc_bndry,dx);
    visc_op->maxOrder(max_order);

    int usehoop = comp == Xvel && (CoordSys::IsRZ());
    int useden  = rho_flag == 1;
    //
    // alpha should be the same size as volume.
    //
    MultiFab alpha(grids,1,GEOM_GROW);

    for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
    {
        DependentMultiFabIterator volumemfi(alphamfi, volume);
        DependentMultiFabIterator rho_halfmfi(alphamfi, (*rho_half));
        assert(alpha.box(alphamfi.index()) == alphamfi.validbox());
        assert(volume.box(alphamfi.index()) == volumemfi.validbox());

        const Box& bx = alphamfi.validbox();

        Array<Real> rcen(bx.length(0));
        parent->Geom(level).GetCellLoc(rcen, bx, 0);

        const int*  lo      = bx.loVect();
        const int*  hi      = bx.hiVect();
        Real*       dat     = alphamfi().dataPtr();
        Box         abx     = ::grow(bx,alpha.nGrow());
        const int*  alo     = abx.loVect();
        const int*  ahi     = abx.hiVect();
        const Real* rcendat = rcen.dataPtr();
        const Real* voli    = volumemfi().dataPtr();
        Box         vbox    = ::grow(volumemfi.validbox(),volume.nGrow());
        const int*  vlo     = vbox.loVect();
        const int*  vhi     = vbox.hiVect();

        const FArrayBox& Rh = rho_halfmfi();
        DEF_CLIMITS(Rh,rho_dat,rlo,rhi);

        FORT_SETALPHA(dat, ARLIM(alo), ARLIM(ahi),
                      lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                      voli, ARLIM(vlo), ARLIM(vhi),
                      rho_dat,ARLIM(rlo),ARLIM(rhi),&usehoop,&useden);
    }

    if (rho_flag == 2)
    {
        MultiFab& S = caller->get_data(State_Type,time);

        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator Smfi(alphamfi, S);
            assert(grids[alphamfi.index()] == alphamfi.validbox());
            alphamfi().mult(Smfi(),alphamfi.validbox(),Density,0,1);
        }
    }
    if (alpha_in != 0)
    {
        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator alpha_inmfi(alphamfi, (*alpha_in));
            assert(grids[alphamfi.index()] == alphamfi.validbox());
            alphamfi().mult(alpha_inmfi(),alphamfi.validbox(),dataComp,0,1);
        }
    }
    if (rhsscale != 0)
    {
        *rhsscale = scale_abec ? 1.0/alpha.max(0) : 1.0;

        visc_op->setScalars(a*(*rhsscale),b*(*rhsscale));
    }
    else
    {
        visc_op->setScalars(a,b);
    }
    if (allnull)
    {
        visc_op->aCoefficients(alpha);
        for (int n = 0; n < BL_SPACEDIM; n++)
        {
            MultiFab bcoeffs(area[n].boxArray(),1,0);
            bcoeffs.copy(area[n]);
            bcoeffs.mult(dx[n]);
            visc_op->bCoefficients(bcoeffs,n);
        }
    }
    else
    {
        visc_op->aCoefficients(alpha);
        for (int n = 0; n < BL_SPACEDIM; n++)
        {
            MultiFab bcoeffs(area[n].boxArray(),1,0);
            bcoeffs.copy(area[n]);
            for (MultiFabIterator bcoeffsmfi(bcoeffs); bcoeffsmfi.isValid();
                 ++bcoeffsmfi)
            {
                DependentMultiFabIterator betanmfi(bcoeffsmfi, (*beta[n]));
                bcoeffsmfi().mult(betanmfi(),dataComp,0,1);
                bcoeffsmfi().mult(dx[n]);
            }
            visc_op->bCoefficients(bcoeffs,n);
        }
    }

    return visc_op;
}

ABecLaplacian*
Diffusion::getViscOp (int                    comp,
                      Real                   a,
                      Real                   b,
                      const MultiFab*        rho_half,
                      int                    rho_flag,
                      Real*                  rhsscale,
                      int                    dataComp,
                      const MultiFab* const* beta,
                      const MultiFab*        alpha_in)
{
    //
    // Note: This assumes that the "NEW" density is to be used, if rho_flag==2
    //
    int allnull, allthere;
    checkBeta(beta, allthere, allnull);

    const Real*  dx     = caller->Geom().CellSize();
    const Box&   domain = caller->Geom().Domain();
    const BCRec& bc     = caller->get_desc_lst()[State_Type].getBC(comp);

    IntVect ref_ratio = level > 0 ? parent->refRatio(level-1) : IntVect::TheUnitVector();

    ViscBndry bndry(grids,1,domain);
    bndry.setHomogValues(bc, ref_ratio);

    ABecLaplacian* visc_op = new ABecLaplacian(bndry,dx);
    visc_op->maxOrder(max_order);

    int usehoop = ((comp==Xvel) && (CoordSys::IsRZ()));
    int useden  = (rho_flag == 1);
    //
    // alpha should be the same size as volume.
    //
    MultiFab alpha(grids,1,GEOM_GROW);

    for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
    {
        DependentMultiFabIterator volumemfi(alphamfi, volume);
        DependentMultiFabIterator rho_halfmfi(alphamfi, (*rho_half));
        assert(alpha.box(alphamfi.index()) == alphamfi.validbox());
        assert(volume.box(alphamfi.index()) == volumemfi.validbox());

        const Box&       bx      = alphamfi.validbox();
        Array<Real> rcen(bx.length(0));
        parent->Geom(level).GetCellLoc(rcen, bx, 0);
        const int*       lo      = bx.loVect();
        const int*       hi      = bx.hiVect();
        Real*            dat     = alphamfi().dataPtr();
        Box              abx     = ::grow(bx,alpha.nGrow());
        const int*       alo     = abx.loVect();
        const int*       ahi     = abx.hiVect();
        const Real*      rcendat = rcen.dataPtr();
        const Real*      voli    = volumemfi().dataPtr();
        Box              vbox    = ::grow(volumemfi.validbox(),volume.nGrow());
        const int*       vlo     = vbox.loVect();
        const int*       vhi     = vbox.hiVect();
        const FArrayBox& Rh      = rho_halfmfi();
        DEF_CLIMITS(Rh,rho_dat,rlo,rhi);

        FORT_SETALPHA(dat, ARLIM(alo), ARLIM(ahi),
                      lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                      voli, ARLIM(vlo), ARLIM(vhi),
                      rho_dat,ARLIM(rlo),ARLIM(rhi),&usehoop,&useden);
    }

    if (rho_flag == 2)
    {
        MultiFab& S = caller->get_new_data(State_Type);

        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator Smfi(alphamfi, S);
            assert(grids[alphamfi.index()] == alphamfi.validbox());
            alphamfi().mult(Smfi(),alphamfi.validbox(),Density,0,1);
        }
    }
    if (alpha_in != 0)
    {
        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator alpha_inmfi(alphamfi, (*alpha_in));
            assert(grids[alphamfi.index()] == alphamfi.validbox());
            alphamfi().mult(alpha_inmfi(),alphamfi.validbox(),dataComp,0,1);
        }
    }
    if (rhsscale != 0)
    {
        *rhsscale = scale_abec ? 1.0/alpha.max(0) : 1.0;

        visc_op->setScalars(a*(*rhsscale),b*(*rhsscale));
    }
    else
    {
        visc_op->setScalars(a,b);
    }
    visc_op->aCoefficients(alpha);

    if (allnull)
    {
        for (int n = 0; n < BL_SPACEDIM; n++)
        {
            MultiFab bcoeffs(area[n].boxArray(),1,0);
            bcoeffs.copy(area[n]);
            bcoeffs.mult(dx[n],0,1,0);
            visc_op->bCoefficients(bcoeffs,n);
        }
    }
    else
    {
        for (int n = 0; n < BL_SPACEDIM; n++)
        {
            MultiFab bcoeffs(area[n].boxArray(),1,0);
            bcoeffs.copy(area[n]);
            for (MultiFabIterator bcoeffsmfi(bcoeffs); bcoeffsmfi.isValid();
                 ++bcoeffsmfi)
            {
                DependentMultiFabIterator betanmfi(bcoeffsmfi, (*beta[n]));
                bcoeffsmfi().mult(betanmfi(),dataComp,0,1);
            }
            bcoeffs.mult(dx[n],0,1,0);
            visc_op->bCoefficients(bcoeffs,n);
        }
    }

    return visc_op;
}

void
Diffusion::getViscTerms (MultiFab&              visc_terms,
                         int                    src_comp,
                         int                    comp, 
                         Real                   time,
                         int                    rho_flag,
                         int                    dataComp,
                         const MultiFab* const* beta)
{
    int allnull, allthere;
    checkBeta(beta, allthere, allnull);
    //
    // Before computing the godunov predictors we may have to
    // precompute the viscous source terms.  To do this we must
    // construct a Laplacian operator, set the coeficients and apply
    // it to the time N data.  First, however, we must precompute the
    // fine N bndry values.  We will do this for each scalar that diffuses.
    //
    // Note: This routine DOES NOT fill grow cells
    //
    const Real* dx = caller->Geom().CellSize();
    MultiFab&   S  = caller->get_data(State_Type,time);

    visc_terms.setVal(0.0,comp-src_comp,1,1);
    //
    // FIXME
    // LinOp classes cannot handle multcomponent MultiFabs yet,
    // construct the components one at a time and copy to visc_terms.
    //
    MultiFab visc_tmp(grids,1,1);
    MultiFab s_tmp(grids,1,1);

    if (is_diffusive[comp])
    {
        ViscBndry visc_bndry;
        getBndryData(visc_bndry,comp,1,time,rho_flag);
        //
        // Set up operator and apply to compute viscous terms.
        //
        const Real a = 0.0;
        const Real b = allnull ? -visc_coef[comp] : -1.0;

        ABecLaplacian visc_op(visc_bndry,dx);

        visc_op.setScalars(a,b);
        visc_op.maxOrder(max_order);

        if (allnull)
        {
            for (int n = 0; n < BL_SPACEDIM; n++)
            {
                MultiFab bcoeffs(area[n].boxArray(),1,0);
                bcoeffs.copy(area[n]);
                bcoeffs.mult(dx[n]);
                visc_op.bCoefficients(bcoeffs,n);
            }
        }
        else
        {
            for (int n = 0; n < BL_SPACEDIM; n++)
            {
                MultiFab bcoeffs(area[n].boxArray(),1,0);
                bcoeffs.copy(area[n]);
                for (MultiFabIterator bcoeffsmfi(bcoeffs);
                     bcoeffsmfi.isValid(); ++bcoeffsmfi)
                {
                    DependentMultiFabIterator betanmfi(bcoeffsmfi, (*beta[n]));
                    bcoeffsmfi().mult(betanmfi(),dataComp,0,1);
                    bcoeffsmfi().mult(dx[n]);
                }
                visc_op.bCoefficients(bcoeffs,n);
            }
        }
        //
        // Copy to single component multifab for operator classes.
        //
        MultiFab::Copy(s_tmp,S,comp,0,1,0);

        if (rho_flag == 2)
        {
            //
            // We want to evaluate (div beta grad) S, not rho*S.
            //
            for (MultiFabIterator Smfi(S); Smfi.isValid(); ++Smfi)
            {
                DependentMultiFabIterator s_tmpmfi(Smfi, s_tmp);
                const Box& box = Smfi.validbox();
                assert(Smfi().min(box,Density) > 0.0);
                s_tmpmfi().divide(Smfi(),box,Density,0,1);
            }
        }

        visc_op.apply(visc_tmp,s_tmp);
        //
        // Must divide by volume.
        //
        for (MultiFabIterator visc_tmpmfi(visc_tmp);
             visc_tmpmfi.isValid(); ++visc_tmpmfi)
        {
            DependentMultiFabIterator volumemfi(visc_tmpmfi, volume);
            assert(grids[visc_tmpmfi.index()] == visc_tmpmfi.validbox());
            visc_tmpmfi().divide(volumemfi(),visc_tmpmfi.validbox(),0,0,1);
        }

#if (BL_SPACEDIM == 2)
        if (comp == Xvel && CoordSys::IsRZ())
        {
            for (MultiFabIterator visc_tmpmfi(visc_tmp);
                 visc_tmpmfi.isValid(); ++visc_tmpmfi)
            {
                DependentMultiFabIterator s_tmpmfi(visc_tmpmfi, s_tmp);
                assert(visc_tmp.box(visc_tmpmfi.index()) == visc_tmpmfi.validbox());
                assert(s_tmp.box(visc_tmpmfi.index()) == s_tmpmfi.validbox());
                //visc_tmp[k] += -mu * u / r^2
                const Box& bx  = visc_tmpmfi.validbox();
                Box        vbx = ::grow(bx,visc_tmp.nGrow());
                Box        sbx = ::grow(s_tmpmfi.validbox(),s_tmp.nGrow());
                Array<Real> rcen(bx.length(0));
                parent->Geom(level).GetCellLoc(rcen, bx, 0);
                const int*  lo      = bx.loVect();
                const int*  hi      = bx.hiVect();
                const int*  vlo     = vbx.loVect();
                const int*  vhi     = vbx.hiVect();
                const int*  slo     = sbx.loVect();
                const int*  shi     = sbx.hiVect();
                Real*       vdat    = visc_tmpmfi().dataPtr();
                Real*       sdat    = s_tmpmfi().dataPtr();
                const Real* rcendat = rcen.dataPtr();
                const Real  mu      = visc_coef[comp];
                FORT_HOOPSRC(ARLIM(lo), ARLIM(hi),
                             vdat, ARLIM(vlo), ARLIM(vhi),
                             sdat, ARLIM(slo), ARLIM(shi),
                             rcendat, &mu);
            }
        }
#endif

        MultiFab::Copy(visc_terms,visc_tmp,0,comp-src_comp,1,0);
    }
}

void
Diffusion::getTensorViscTerms (MultiFab&              visc_terms, 
                               Real                   time,
                               int                    dataComp,
                               const MultiFab* const* beta)
{
    int allthere;
    checkBeta(beta, allthere);

    const int src_comp = Xvel;
    const int ncomp    = visc_terms.nComp();

    if (ncomp < BL_SPACEDIM)
        BoxLib::Abort("Diffusion::getTensorViscTerms(): visc_terms needs at least BL_SPACEDIM components");

    //
    // Before computing the godunov predicitors we may have to
    // precompute the viscous source terms.  To do this we must
    // construct a Laplacian operator, set the coeficients and apply
    // it to the time N data.  First, however, we must precompute the
    // fine N bndry values.  We will do this for each scalar that diffuses.
    //
    // Note: This routine DOES NOT fill grow cells
    //
    const Real* dx   = caller->Geom().CellSize();
    MultiFab&   S    = caller->get_data(State_Type,time);
    const int   ngrd = grids.length();
    visc_terms.setVal(0.0,src_comp,BL_SPACEDIM,1);
    //
    // FIXME
    // LinOp classes cannot handle multcomponent MultiFabs yet,
    // construct the components one at a time and copy to visc_terms.
    //
    MultiFab visc_tmp(grids,BL_SPACEDIM,1);
    MultiFab s_tmp(grids,BL_SPACEDIM,1);

    if (is_diffusive[src_comp])
    {
        ViscBndryTensor visc_bndry;
        getTensorBndryData(visc_bndry,time);
        //
        // Set up operator and apply to compute viscous terms.
        //
        const Real a = 0.0;
        const Real b = -1.0;

        DivVis tensor_op(visc_bndry,dx);
        tensor_op.maxOrder(tensor_max_order);
        tensor_op.setScalars(a,b);

        const int nghost = 0;
        //
        // alpha should be the same size as volume.
        //
        const int nCompAlpha = BL_SPACEDIM == 2  ?  2 : 1;
        MultiFab alpha(grids,nCompAlpha,nghost);
        alpha.setVal(0.0);
        tensor_op.aCoefficients(alpha);

        for (int n = 0; n < BL_SPACEDIM; n++)
        {
            MultiFab bcoeffs(area[n].boxArray(),1,nghost);
            bcoeffs.setBndry(0);
            MultiFab::Copy(bcoeffs,area[n],0,0,1,nghost);
            for (MultiFabIterator bcoeffsmfi(bcoeffs); bcoeffsmfi.isValid();
                 ++bcoeffsmfi)
            {
                DependentMultiFabIterator betanmfi(bcoeffsmfi, (*beta[n]));
                bcoeffsmfi().mult(dx[n]);
                bcoeffsmfi().mult(betanmfi(),dataComp,0,1);
            }
            tensor_op.bCoefficients(bcoeffs,n);
        }

        MultiFab::Copy(s_tmp,S,Xvel,0,BL_SPACEDIM,0);

        tensor_op.apply(visc_tmp,s_tmp);
        //
        // Must divide by volume.
        //
        for (MultiFabIterator visc_tmpmfi(visc_tmp);
             visc_tmpmfi.isValid(); ++visc_tmpmfi)
        {
            DependentMultiFabIterator volumemfi(visc_tmpmfi, volume);
            assert(grids[volumemfi.index()] == volumemfi.validbox());
            visc_tmpmfi().divide(volumemfi(),volumemfi.validbox(),0,0,1);
        }
#if (BL_SPACEDIM == 2)
        if (CoordSys::IsRZ())
        {
            int fort_xvel_comp = Xvel+1;
            for (int k = 0; k < ngrd; k++)
            {
                //visc_tmp[k] += -mu * u / r^2
                const Box& bx  = visc_tmp.box(k);
                Box        vbx = ::grow(bx,visc_tmp.nGrow());
                Box        sbx = ::grow(s_tmp.box(k),s_tmp.nGrow());
                Array<Real> rcen(bx.length(0));
                parent->Geom(level).GetCellLoc(rcen, bx, 0);
                const int*       lo        = bx.loVect();
                const int*       hi        = bx.hiVect();
                const int*       vlo       = vbx.loVect();
                const int*       vhi       = vbx.hiVect();
                const int*       slo       = sbx.loVect();
                const int*       shi       = sbx.hiVect();
                Real*            vdat      = visc_tmp[k].dataPtr();
                Real*            sdat      = s_tmp[k].dataPtr();
                const Real*      rcendat   = rcen.dataPtr();
                const FArrayBox& betax     = (*beta[0])[k];
                const Real*      betax_dat = betax.dataPtr(dataComp);
                const int*       betax_lo  = betax.loVect();
                const int*       betax_hi  = betax.hiVect();
                const FArrayBox& betay     = (*beta[1])[k];
                const Real*      betay_dat = betay.dataPtr(dataComp);
                const int*       betay_lo  = betay.loVect();
                const int*       betay_hi  = betay.hiVect();

                FORT_TENSOR_HOOPSRC(&fort_xvel_comp,ARLIM(lo), ARLIM(hi),
                                    vdat, ARLIM(vlo), ARLIM(vhi),
                                    sdat, ARLIM(slo), ARLIM(shi),
                                    rcendat, 
                                    betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                                    betay_dat,ARLIM(betay_lo),ARLIM(betay_hi));
            }
        }
#endif

        MultiFab::Copy(visc_terms,visc_tmp,0,0,BL_SPACEDIM,0);
    }
}

void
Diffusion::getBndryData (ViscBndry& bndry,
                         int        src_comp,
                         int        num_comp,
                         Real       time,
                         int        rho_flag)
{
    assert(num_comp == 1);
    //
    // Fill phys bndry vals of grow cells of (tmp) multifab passed to bndry.
    //
    // TODO -- A MultiFab is a huge amount of space in which to pass along
    // the phys bc's.  InterpBndryData needs a more efficient interface.
    //
    //
    const int    nGrow = 1;
    const BCRec& bc    = caller->get_desc_lst()[State_Type].getBC(src_comp);

    MultiFab S(grids, num_comp, nGrow);

    bndry.define(grids,num_comp,caller->Geom());

    FillPatchIterator Phi_fpi(*caller,S,nGrow,time,State_Type,src_comp,num_comp);

    FillPatchIterator Rho_fpi(*caller,S);

    if (rho_flag == 2)
	Rho_fpi.Initialize(nGrow,time,State_Type,Density,1);

    for ( ; Phi_fpi.isValid(); ++Phi_fpi)
    {
	if (rho_flag == 2)
	    Rho_fpi.isValid();

        DependentMultiFabIterator S_mfi(Phi_fpi, S);

        const BoxList gCells = ::boxDiff(Phi_fpi().box(), Phi_fpi.validbox());

        for (BoxListIterator bli(gCells); bli; ++bli)
        {
            S_mfi().copy(Phi_fpi(),bli(),0,bli(),0,num_comp);

            if (rho_flag == 2)
                for (int n = 0; n < num_comp; ++n)
                    S_mfi().divide(Rho_fpi(),bli(),0,n,1);
        }

	if (rho_flag == 2)
	    ++Rho_fpi;
    }
    
    if (level == 0)
    {
        bndry.setBndryValues(S,0,0,num_comp,bc);
    }
    else
    {
        BoxArray cgrids = grids;
        cgrids.coarsen(crse_ratio);
        BndryRegister crse_br(cgrids,0,1,1,num_comp);
        //
        // interp for solvers over ALL c-f brs, need safe data.
        //
        crse_br.setVal(BL_BOGUS);
        coarser->FillBoundary(crse_br,src_comp,0,num_comp,time,rho_flag);
        bndry.setBndryValues(crse_br,0,S,0,0,num_comp,crse_ratio,bc);
    }
}

void
Diffusion::FillBoundary (BndryRegister& bdry,
                         int            src_comp,
                         int            dest_comp,
                         int            num_comp,
                         Real           time,
                         int            rho_flag)
{
    //
    // Need one grow cell filled for linear solvers.
    // We assume filPatch gets this right, where possible.
    //
    const int nGrow = 1;

    MultiFab S(caller->boxArray(),num_comp,nGrow);

    {
        FillPatchIterator S_fpi(*caller,S,nGrow,time,State_Type,src_comp,num_comp);

        for ( ; S_fpi.isValid(); ++S_fpi)
        {
            S[S_fpi.index()].copy(S_fpi(), 0, 0, num_comp);
        }
    }

    if (rho_flag == 2)
    {
        FillPatchIterator rho_fpi(*caller,S,nGrow,time,State_Type,Density,1);

        for ( ; rho_fpi.isValid(); ++rho_fpi)
        {
            S[rho_fpi.index()].divide(rho_fpi());
        }
    }
    //
    // Copy into boundary register.
    //
    bdry.copyFrom(S,nGrow,0,dest_comp,num_comp);
}

void
Diffusion::getTensorBndryData (ViscBndryTensor& bndry, 
                               Real             time)
{
    const int num_comp = BL_SPACEDIM;
    const int src_comp = Xvel;
    const int nDer     = MCLinOp::bcComponentsNeeded();
    //
    // Create the BCRec's interpreted by ViscBndry objects
    //
    Array<BCRec> bcarray(nDer, BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
                                     D_DECL(EXT_DIR,EXT_DIR,EXT_DIR)));

    for (int idim = 0; idim < BL_SPACEDIM; idim++)
        bcarray[idim] = caller->get_desc_lst()[State_Type].getBC(src_comp+idim);

    bndry.define(grids,nDer,caller->Geom());

    const int nGrow = 1;
    MultiFab S(grids,num_comp,nGrow,Fab_allocate);
    FillPatchIterator Phi_fpi(*caller,S,nGrow,time,State_Type,src_comp,num_comp);
    for ( ; Phi_fpi.isValid(); ++Phi_fpi)
    {
        DependentMultiFabIterator S_mfi(Phi_fpi, S);
        
        const BoxList gCells = ::boxDiff(Phi_fpi().box(), Phi_fpi.validbox());
        
        for (BoxListIterator bli(gCells); bli; ++bli)
        {
            S_mfi().copy(Phi_fpi(),bli(),0,bli(),0,num_comp);
        }
    }
    
    if (level == 0)
    {
        bndry.setBndryValues(S,0,0,num_comp,bcarray);
    }
    else
    {
        BoxArray cgrids(grids);
        cgrids.coarsen(crse_ratio);
        BndryRegister crse_br(cgrids,0,1,1,num_comp);
        crse_br.setVal(BL_BOGUS);
        const int rho_flag = 0;
        coarser->FillBoundary(crse_br,src_comp,0,num_comp,time,rho_flag);
        bndry.setBndryValues(crse_br,0,S,0,0,num_comp,crse_ratio[0],bcarray);
    }
}

void
Diffusion::checkBetas (const MultiFab* const* beta1, 
                       const MultiFab* const* beta2,
                       int&                   allthere,
                       int&                   allnull) const
{
    int allnull1, allnull2, allthere1, allthere2;

    checkBeta(beta1,allthere1,allnull1);
    checkBeta(beta2,allthere2,allnull2);
    allnull  = allnull1 && allnull2;
    allthere = allthere1 && allthere2;

    if (!(allthere || allnull))
        BoxLib::Abort("Diffusion::checkBetas(): betas must either be all 0 or all non-0");
}

void
Diffusion::checkBeta (const MultiFab* const* beta,
                      int&                   allthere,
                      int&                   allnull) const
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
        BoxLib::Abort("Diffusion::checkBeta(): betas must be all 0 or all non-0");
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
        BoxLib::Abort("Diffusion::checkBeta(): betas must be all non-0");
}

void
Diffusion::allocFluxBoxesLevel (MultiFab**& fluxbox, 
                                int         nghost,
                                int         nvar)
{
    fluxbox = new MultiFab*[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        BoxArray edge_boxes(grids);
        edge_boxes.surroundingNodes(dir);
        fluxbox[dir] = new MultiFab(edge_boxes,nvar,nghost);
    }
}

void
Diffusion::removeFluxBoxesLevel (MultiFab**& fluxbox) 
{
    if (fluxbox != 0)
    {
        for (int i = 0; i<BL_SPACEDIM; i++)
            delete fluxbox[i];
        delete [] fluxbox;
        fluxbox = 0;
    }
}

#ifdef USE_NAVIERSTOKES
//
// This routine computes the vector div mu SI, where I is the identity 
// tensor, S = div U, and mu is constant.
//

void
Diffusion::compute_divmusi (Real      time,
                            Real      mu,
                            MultiFab& divmusi)
{
    //
    // Compute on valid region.  Calling function should fill grow cells.
    //
    if (mu > 0.0)
    {
        const int     nGrowDU  = 1;
        const Real*   dx       = caller->Geom().CellSize();
        NavierStokes& ns_level = *(NavierStokes*) &(parent->getLevel(level));
        MultiFab*     divu_fp  = ns_level.getDivCond(nGrowDU,time);

        for (MultiFabIterator divmusimfi(divmusi); divmusimfi.isValid(); ++divmusimfi)
        {
            DependentMultiFabIterator divumfi(divmusimfi,*divu_fp);

            FArrayBox& fab  = divmusimfi();
            FArrayBox& divu = divumfi();
            const Box& box  = divmusimfi.validbox();

            FORT_DIV_MU_SI(box.loVect(), box.hiVect(), dx, &mu,
                           ARLIM(divu.loVect()), ARLIM(divu.hiVect()),
                           divu.dataPtr(),
                           ARLIM(fab.loVect()),  ARLIM(fab.hiVect()),
                           fab.dataPtr());
        }

        delete divu_fp;
    }
    else
    {
        const int nGrow = 0; // Not to fill grow cells here
        divmusi.setVal(0.0,nGrow);
    }
}

//
// This routine computes the vector div beta SI, where I is the identity 
// tensor, S = div U, and beta is non-constant.
//

void
Diffusion::compute_divmusi (Real                   time,
                            const MultiFab* const* beta,
                            MultiFab&              divmusi)
{
    const int     nGrowDU  = 1;
    const Real*   dx       = caller->Geom().CellSize();
    NavierStokes& ns_level = *(NavierStokes*) &(parent->getLevel(level));
    MultiFab*     divu_fp  = ns_level.getDivCond(nGrowDU,time);

    for (MultiFabIterator divmusimfi(divmusi); divmusimfi.isValid(); ++divmusimfi)
    {
        DependentMultiFabIterator divumfi(divmusimfi,*divu_fp);

        const int  i    = divmusimfi.index();
        FArrayBox& fab  = divmusimfi();
        FArrayBox& divu = divumfi();
        const Box& box  = divmusimfi.validbox();

        DependentMultiFabIterator beta0mfi(divmusimfi, (*beta[0]));
        DependentMultiFabIterator beta1mfi(divmusimfi, (*beta[1]));

        DEF_CLIMITS(beta0mfi(),betax,betaxlo,betaxhi);
        DEF_CLIMITS(beta1mfi(),betay,betaylo,betayhi);

#if (BL_SPACEDIM==3)
        DependentMultiFabIterator beta2mfi(divmusimfi, (*beta[2]));
        DEF_CLIMITS(beta2mfi(),betaz,betazlo,betazhi);
#endif
        assert(grids[divmusimfi.index()] == divmusimfi.validbox());

        FORT_DIV_VARMU_SI(box.loVect(),box.hiVect(), dx,
                          ARLIM(divu.loVect()), ARLIM(divu.hiVect()),
                          divu.dataPtr(),
                          ARLIM(betaxlo), ARLIM(betaxhi), betax,
                          ARLIM(betaylo), ARLIM(betayhi), betay,
#if (BL_SPACEDIM==3)
                          ARLIM(betazlo), ARLIM(betazhi), betaz,
#endif
                          ARLIM(divmusi[i].loVect()), ARLIM(divmusi[i].hiVect()),
                          divmusi[i].dataPtr());
    }

    delete divu_fp;
}
#endif /*USE_NAVIERSTOKES*/
