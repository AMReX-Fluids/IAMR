//BL_COPYRIGHT_NOTICE

//
// $Id: Diffusion.cpp,v 1.46 1998-11-19 23:38:08 lijewski Exp $
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
const Real BL_SAFE_BOGUS = -666.e30;
#else
const REAL BL_BOGUS      = FLT_QNAN;
const Real BL_SAFE_BOGUS = -666.e30;
#endif
#else
const REAL BL_BOGUS      = 1.e30;
const Real BL_SAFE_BOGUS = -666.e30;
#endif

#if (BL_SPACEDIM==2) && defined (USE_TENSOR)
//
// Include files for tensor solve.
//
#include <DivVis.H>
#include <LO_BCTYPES.H>
#include <MCMultiGrid.H>
#include <MCCGSolver.H>
#include <ViscBndry2D.H>
#endif

const char NL = '\n';

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

int  Diffusion::first = 1;
int  Diffusion::do_reflux = 1;
int  Diffusion::use_cg_solve = 0;
int  Diffusion::use_tensor_cg_solve = 0;
bool Diffusion::use_mg_precond_flag = false;
int  Diffusion::verbose = 0;
int  Diffusion::max_order = 2;
int  Diffusion::tensor_max_order = 2;
int  Diffusion::use_dv_constant_mu_def = 1;
int  Diffusion::scale_abec = 0;
int  Diffusion::est_visc_mag = 1;
int  Diffusion::Lphi_in_abs_tol = 0;

Array<Real> Diffusion::typical_vals;
const Real typical_vals_DEF = 1.0;

//
// This modules own MultiFab - MultiFab copy with ghost cells
// between MultiFabs with the same ProcessorMap().
//

static
void
Copy (MultiFab&       dst,
      const MultiFab& src,
      int             srccomp,
      int             dstcomp,
      int             numcomp,
      int             nghost)
{
    assert(dst.nGrow() >= nghost && src.nGrow() >= nghost);

    for (MultiFabIterator mfi(dst); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi,src);

        Box bx = ::grow(mfi.validbox(),nghost) & ::grow(dmfi.validbox(),nghost);

        if (bx.ok())
        {
            mfi().copy(dmfi(),
                       bx,
                       srccomp,
                       bx,
                       dstcomp,
                       numcomp);
        }
    }
}

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
        ppdiff.query("use_dv_constant_mu",use_dv_constant_mu_def);
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

        int n_visc = _visc_coef.length();
        int n_diff = _is_diffusive.length();
        if (n_diff < NUM_STATE || n_visc < NUM_STATE)
        {
            cout << "Diffusion::Diffusion(): is_diffusive and/or visc_coef arrays are " <<
                " not long enough\n";
            BoxLib::Abort("Bye");
        }
        visc_coef.resize(NUM_STATE);
        is_diffusive.resize(NUM_STATE);
        for (int i = 0; i < NUM_STATE; i++) {
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

Diffusion::~Diffusion() {}

void
Diffusion::echo_settings () const
{
    //
    // Print out my settings.
    //
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "Diffusion settings..." << NL;
        cout << "  From diffuse:" << NL;
        cout << "   use_cg_solve =        " << use_cg_solve << NL;
        cout << "   use_tensor_cg_solve = "  << use_tensor_cg_solve << NL;
        cout << "   use_dv_constant_mu =  " << use_dv_constant_mu_def << NL;
        cout << "   use_mg_precond_flag = " << use_mg_precond_flag << NL;
        cout << "   max_order =           " << max_order << NL;
        cout << "   tensor_max_order =    " << tensor_max_order << NL;
        cout << "   scale_abec =          " << scale_abec << NL;
        cout << "   est_visc_mag =        " << est_visc_mag << NL;
        cout << "   Lphi_in_abs_tol =     " << Lphi_in_abs_tol << NL;
    
        cout << "   typical_vals =";
        for (int i = 0; i <NUM_STATE; i++)
        {
            cout << "  " << typical_vals[i];
        }
        cout << NL;
    
        cout << NL;
        cout << "  From ns:" << NL;
        cout << "   do_reflux =           " << do_reflux << NL;
        cout << "   visc_tol =            " << visc_tol << NL;
        cout << "   visc_abs_tol =        " << visc_abs_tol << NL;
    
        cout << "   is_diffusive =";
        for (int i =0; i < NUM_STATE; i++)
        {
            cout << "  " << is_diffusive[i];
        }
        cout << NL;
    
        cout << "   visc_coef =";
        for (int i = 0; i < NUM_STATE; i++)
        {
            cout << "  " << visc_coef[i];
        }
        cout << NL;
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
    
	int do_const_visc = (!est_visc_mag) || 
                            ( (!allthere_n) && (!allthere_np1) ); 
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
Diffusion::diffuse_scalar (Real       dt,
                           int        sigma,
                           Real       be_cn_theta,
                           MultiFab*  rho_half,
                           int        rho_flag,
                           int        do_viscreflux,
                           MultiFab*  delta_rhs, 
                           MultiFab*  alpha, 
                           MultiFab** betan, 
                           MultiFab** betanp1)
{
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "... diffuse_scalar " << sigma << NL;
    }

    int allnull, allthere;
    checkBetas(betan, betanp1, allthere, allnull);

    const int finest_level = parent->finestLevel();
    //
    // At this point, S_old has bndry at time N.
    //
    MultiFab& S_old = caller->get_old_data(State_Type);
    MultiFab& S_new = caller->get_new_data(State_Type);

    MultiFab* Rho_old = 0;
    MultiFab* Rho_new = 0;

    if (rho_flag == 2)
    {
        Rho_old = new MultiFab(grids,1,1,Fab_allocate);
        Copy(*Rho_old,S_old,Density,0,1,1);
        Rho_new = new MultiFab(grids,1,1,Fab_allocate);
        Copy(*Rho_new,S_new,Density,0,1,1);
    }

    const Box& domain    = caller->Geom().Domain();
    const int* domlo     = domain.loVect();
    const int* domhi     = domain.hiVect();
    const Real* dx       = caller->Geom().CellSize();
    const Real cur_time  = caller->get_state_data(State_Type).curTime();
    const Real prev_time = caller->get_state_data(State_Type).prevTime();

    FArrayBox xflux,yflux,zflux;
    //
    // We are now at this point in NS::scalar_update:
    //if (is_diffusive[sigma]) {
    //
    // S_new now contains s(n) - dt*aofs.  We will use
    // this as part of the RHS for the viscous solve.
    //
    MultiFab Rhs(grids,1,0,Fab_allocate);
    MultiFab Soln(grids,1,1,Fab_allocate);
    {
        //
        // Set up Rhs.
        //
        const Real a = 0.0;
        Real b       = -(1.0-be_cn_theta)*dt;
        if (allnull)
            b *= visc_coef[sigma];

        ViscBndry visc_bndry;
        ABecLaplacian* visc_op = getViscOp(sigma,a,b,prev_time,visc_bndry,rho_half,rho_flag,0,betan);
        visc_op->maxOrder(max_order);

        if (rho_flag == 2)
        {
            //
            // We are going to solve for S, not rho*S.
            //
            for (MultiFabIterator S_oldmfi(S_old); S_oldmfi.isValid(); ++S_oldmfi)
            {
                DependentMultiFabIterator Rho_oldmfi(S_oldmfi, (*Rho_old));
                assert(grids[S_oldmfi.index()] == S_oldmfi.validbox());
                Box box = ::grow(S_oldmfi.validbox(),1);
                S_oldmfi().divide(Rho_oldmfi(),box,0,sigma,1);
            }
        }
        //
        // Copy to single-component multifab.
        // Note: use Soln as a temporary here.
        //
        Copy(Soln,S_old,sigma,0,1,1);

        visc_op->apply(Rhs,Soln);

        delete visc_op;
        //
        // Copy back to S_old for use in creating viscous fluxes.
        // NOTE: this requires that the visc_op->apply call returns
        //       Soln with the ghost cells correctly filled.
        //
        Copy(S_old,Soln,0,sigma,1,1);
        //
        // Complete Rhs by adding body sources.
        //
        for (MultiFabIterator S_newmfi(S_new); S_newmfi.isValid(); ++S_newmfi)
        {
            DependentMultiFabIterator volumemfi(S_newmfi, volume);
            DependentMultiFabIterator rho_halfmfi(S_newmfi, (*rho_half));
            DependentMultiFabIterator Rhsmfi(S_newmfi, Rhs);
            assert(grids[S_newmfi.index()] == S_newmfi.validbox());
            //
            // Scale inviscid part by volume.
            //
            S_newmfi().mult(volumemfi(),S_newmfi.validbox(),0,sigma,1);

            if (rho_flag == 1)
            {
                //
                // Multiply by density at time nph.
                //
                S_newmfi().mult(rho_halfmfi(),S_newmfi.validbox(),0,sigma,1);
            }

            if (alpha!=0)
            {
                DependentMultiFabIterator alphamfi(S_newmfi, (*alpha));
                S_newmfi().mult(alphamfi(),S_newmfi.validbox(),0,sigma,1);
            }
            //
            // Add to rhs.
            //
            Rhsmfi().plus(S_newmfi(),S_newmfi.validbox(),sigma,0,1);
        }
        if (delta_rhs != 0)
        {
            for (MultiFabIterator delta_rhsmfi(*delta_rhs);
                 delta_rhsmfi.isValid(); ++delta_rhsmfi)
            {
                DependentMultiFabIterator volumemfi(delta_rhsmfi, volume);
                DependentMultiFabIterator Rhsmfi(delta_rhsmfi, Rhs);
                assert(grids[delta_rhsmfi.index()] == delta_rhsmfi.validbox());
                delta_rhsmfi().mult(dt);
                delta_rhsmfi().mult(volumemfi(),delta_rhsmfi.validbox(),0,0,1);
                Rhsmfi().plus(delta_rhsmfi(),delta_rhsmfi.validbox(),0,0,1);
            }
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
    ABecLaplacian* visc_op = getViscOp(sigma,a,b,cur_time,visc_bndry,rho_half,rho_flag,&rhsscale,betanp1,alpha);
    Rhs.mult(rhsscale,0,1);
    visc_op->maxOrder(max_order);
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

    int visc_op_lev = 0;
    visc_op->applyBC(Soln,visc_op_lev);

    delete visc_op;
    //
    // Copy into state variable at new time, with bcs hopefully.
    //
    Copy(S_new,Soln,0,sigma,1,1);
    //
    // Create diffusive fluxes here.
    //
    if (do_reflux && do_viscreflux)
    {
        FArrayBox flux;

        for (MultiFabIterator S_oldmfi(S_old); S_oldmfi.isValid(); ++S_oldmfi)
        {
            DependentMultiFabIterator S_newmfi(S_oldmfi, S_new);

            int i                  = S_oldmfi.index();
            const Box& grd         = S_oldmfi.validbox();
            const FArrayBox& sold  = S_oldmfi();
            const FArrayBox& snew  = S_newmfi();

            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                DependentMultiFabIterator areadmfi(S_oldmfi, area[d]);

                flux.resize(::surroundingNodes(grd,d), 1);

                const FArrayBox& aread = areadmfi();

                if (allthere)
                {
                    const Real mult = -1.0;

                    DependentMultiFabIterator betaoldmfi(S_oldmfi,*betan[d]);
                    DependentMultiFabIterator betanewmfi(S_oldmfi,*betanp1[d]);

                    const FArrayBox& betaold = betaoldmfi();
                    const FArrayBox& betanew = betanewmfi();

                    FORT_VISCFLUX_VC(sold.dataPtr(sigma),
                                     snew.dataPtr(sigma),
                                     ARLIM(sold.loVect()), ARLIM(sold.hiVect()),
                                     grd.loVect(), grd.hiVect(),
                                     flux.dataPtr(),
                                     ARLIM(flux.loVect()), ARLIM(flux.hiVect()),
                                     aread.dataPtr(),
                                     ARLIM(aread.loVect()),ARLIM(aread.hiVect()),
                                     betaold.dataPtr(), betanew.dataPtr(),
                                     ARLIM(betaold.loVect()),ARLIM(betaold.hiVect()),
                                     &dx[d],&mult,&be_cn_theta,&d);
                }
                else
                {
                    const Real mult = -visc_coef[sigma];

                    FORT_VISCFLUX_CC(sold.dataPtr(sigma),
                                     snew.dataPtr(sigma),
                                     ARLIM(sold.loVect()), ARLIM(sold.hiVect()),
                                     grd.loVect(), grd.hiVect(),
                                     flux.dataPtr(),
                                     ARLIM(flux.loVect()), ARLIM(flux.hiVect()),
                                     aread.dataPtr(),
                                     ARLIM(aread.loVect()),ARLIM(aread.hiVect()),
                                     &dx[d],&mult,&be_cn_theta,&d);
                }
                //
                // Fluxes expected to be in extensive form.
                //
                if (level < finest_level)
                    finer->viscflux_reg->CrseInit(flux,flux.box(),d,0,sigma,1,-dt);

                if (level > 0)
                    viscflux_reg->FineAdd(flux,d,i,0,sigma,1,dt);               
            }
        }
        if (level < finest_level)
            finer->viscflux_reg->CrseInitFinish();
    } 

    if (rho_flag == 2)
    {
        //
        // Return S_old to hold rho*S
        // Also we solved for S, not rho*S, so fix S_new.
        //
        for (MultiFabIterator S_oldmfi(S_old); S_oldmfi.isValid(); ++S_oldmfi)
        {
            DependentMultiFabIterator S_newmfi(S_oldmfi, S_new);
            DependentMultiFabIterator Rho_newmfi(S_oldmfi,*Rho_new);
            DependentMultiFabIterator Rho_oldmfi(S_oldmfi,*Rho_old);

            assert(grids[S_oldmfi.index()] == S_oldmfi.validbox());

            Box box = ::grow(S_oldmfi.validbox(),1);

            S_newmfi().mult(Rho_newmfi(),box,0,sigma,1);
            S_oldmfi().mult(Rho_oldmfi(),box,0,sigma,1);
        }

        delete Rho_old;
        delete Rho_new;
    }
}

void
Diffusion::diffuse_velocity (Real       dt,
                             Real       be_cn_theta,
                             MultiFab*  rho_half,
                             int        rho_flag,
                             MultiFab*  delta_rhs,
                             MultiFab** betan, 
                             MultiFab** betanp1)
{
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "... diffuse_velocity\n";
    }

    int allnull, allthere;
    checkBetas(betan, betanp1, allthere, allnull);

    int constant_viscosity = allnull;
    int use_dv_constant_mu = use_dv_constant_mu_def;

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        use_dv_constant_mu = use_dv_constant_mu && visc_coef[Xvel+i] >= 0.0;
    }

    if (!use_dv_constant_mu && use_dv_constant_mu_def)
    {
        cout << "Diffusion::diffuse_velocity() : must have velocity visc_coefs "
             << ">= 0.0 if use_dv_constant_mu == 1\n";
        BoxLib::Abort("Exiting.");
    }

    if (constant_viscosity || use_dv_constant_mu)
    {
        diffuse_velocity_constant_mu(dt, be_cn_theta, rho_half, delta_rhs);
    }
    else
    {
        diffuse_tensor_velocity(dt,be_cn_theta,rho_half,delta_rhs,betan,betanp1);
    }
}

void
Diffusion::diffuse_velocity_constant_mu (Real      dt,
                                         Real      be_cn_theta,
                                         MultiFab* rho_half,
                                         MultiFab* delta_rhs)
{
    //
    // At this point, S_old has bndry at time N.  S_new contains GRAD(SU).
    //
    MultiFab& U_old        = caller->get_old_data(State_Type);
    MultiFab& U_new        = caller->get_new_data(State_Type);
    const Box& domain      = caller->Geom().Domain();
    const int* domlo       = domain.loVect();
    const int* domhi       = domain.hiVect();
    const Real* dx         = caller->Geom().CellSize();
    const Real cur_time    = caller->get_state_data(State_Type).curTime();
    const Real prev_time   = caller->get_state_data(State_Type).prevTime();
    const int finest_level = parent->finestLevel();

    FArrayBox xflux, yflux, zflux;
    //
    // At this point in time we can only do decoupled scalar
    // so we loop over components.
    //
    for (int comp = 0; comp < BL_SPACEDIM; comp++)
    {
        int sigma = Xvel + comp;
        //
        // U_new now contains the inviscid update of U.
        // This is part of the RHS for the viscous solve.
        //
        MultiFab Rhs(grids,1,0,Fab_allocate);
        MultiFab Soln(grids,1,1,Fab_allocate);

        int rho_flag = 1;
        {
            //
            // Set up Rhs.
            //
            const Real a = 0.0;
            const Real b = -(1.0-be_cn_theta)*visc_coef[sigma]*dt;

            ViscBndry visc_bndry;
            ABecLaplacian* visc_op = getViscOp(sigma,a,b,prev_time,visc_bndry,rho_half,rho_flag);
            visc_op->maxOrder(max_order);
            //
            // Copy to single-component multifab.
            // Note: use Soln as a temporary here.
            //
            Copy(Soln,U_old,sigma,0,1,1);
            visc_op->apply(Rhs,Soln);
            delete visc_op;
            //
            // Complete Rhs by adding body sources.
            //
            for (MultiFabIterator U_newmfi(U_new); U_newmfi.isValid(); ++U_newmfi)
            {
                DependentMultiFabIterator volumemfi(U_newmfi, volume);
                DependentMultiFabIterator Rhsmfi(U_newmfi, Rhs);
                DependentMultiFabIterator rho_halfmfi(U_newmfi, (*rho_half));

                assert(grids[U_newmfi.index()] == U_newmfi.validbox());
                //
                // Scale inviscid part by volume.
                //
                U_newmfi().mult(volumemfi(),U_newmfi.validbox(),0,sigma,1);
                //
                // Multiply by density at time nph.
                //
                U_newmfi().mult(rho_halfmfi(),U_newmfi.validbox(),0,sigma,1);
                //
                // Add to Rhs which contained mu/(2 dt) Lap(u).
                //
                Rhsmfi().plus(U_newmfi(),U_newmfi.validbox(),sigma,0,1);
            }

            if (delta_rhs != 0)
            {
                for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
                {
                    DependentMultiFabIterator volumemfi(Rhsmfi, volume);
                    DependentMultiFabIterator delta_rhsmfi(Rhsmfi,(*delta_rhs));
                    assert(grids[Rhsmfi.index()] == Rhsmfi.validbox());
                    delta_rhsmfi().mult(dt,comp,1);
                    delta_rhsmfi().mult(volumemfi(),Rhsmfi.validbox(),0,comp,1);
                    Rhsmfi().plus(delta_rhsmfi(),Rhsmfi.validbox(),comp,0,1);
                }
            }
            //
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
                    DependentMultiFabIterator U_oldmfi(Rhsmfi, U_old);
                    assert(Rhs.box(Rhsmfi.index()) == Rhsmfi.validbox());
                    assert(U_old.box(Rhsmfi.index()) == U_oldmfi.validbox());
                    assert(volume.box(Rhsmfi.index()) == volumemfi.validbox());

                    const Box& bx = Rhsmfi.validbox();
                    Box sbx       = ::grow(U_oldmfi.validbox(),U_old.nGrow());
                    Array<Real> rcen(bx.length(0));
                    parent->Geom(level).GetCellLoc(rcen, bx, 0);
                    const int* lo       = bx.loVect();
                    const int* hi       = bx.hiVect();
                    const int* slo      = sbx.loVect();
                    const int* shi      = sbx.hiVect();
                    Real* rhs           = Rhsmfi().dataPtr();
                    Real* sdat          = U_oldmfi().dataPtr(sigma);
                    const Real* rcendat = rcen.dataPtr();
                    const Real coeff    = (1.0-be_cn_theta)*visc_coef[sigma]*dt;
                    const Real* voli    = volumemfi().dataPtr();
                    Box vbox            = ::grow(volumemfi.validbox(),volume.nGrow());
                    const int* vlo      = vbox.loVect();
                    const int* vhi      = vbox.hiVect();
                    FORT_HOOPRHS(rhs, ARLIM(lo), ARLIM(hi), 
                                 sdat, ARLIM(slo), ARLIM(shi),
                                 rcendat, &coeff, voli, ARLIM(vlo),ARLIM(vhi));
                }
            }
#endif
        }
        //
        // Compute guess of solution.
        //
        if (level == 0)
        {
            Copy(Soln,U_old,sigma,0,1,1);
        }
        else
        {
            caller->FillCoarsePatch(Soln,0,cur_time,State_Type,sigma,1);
        }
        //
        // Copy guess into U_new.
        //
        // The new-time operator is initialized with a "guess" for the
        // new-time state.  We intentionally initialize the grow cells with
        // a bogus value to emphasize that the values are not to be
        // considered "valid" (we shouldn't specify any grow cell
        // information), but rather are to be filled by the "physics bc's,
        // etc" in the problem-dependent code.  In the course of this
        // filling (typically while generating/filling the BndryData object
        // for the solvers), StateData::filcc is called to get physical
        // bc's.  Here 'something computable' has to already exist in the
        // grow cells (even though filcc ultimately will fill the corner
        // correctly, if applicable).  This is apparently is where the
        // `something computable' is to be set.
        //
        int n_comp  = 1;
        int n_ghost = 1;
        U_new.setVal(BL_SAFE_BOGUS,sigma,n_comp,n_ghost);
        n_ghost = 0;
        U_new.copy(Soln,0,sigma,n_comp);
        //
        // Construct viscous operator with bndry data at time N+1.
        //
        const Real a = 1.0;
        const Real b = be_cn_theta*dt*visc_coef[sigma];
       
        ViscBndry visc_bndry;
        rho_flag      = 1 ;
        Real rhsscale = 1.0;
        ABecLaplacian* visc_op = getViscOp(sigma,a,b,cur_time,visc_bndry,rho_half,rho_flag,&rhsscale);
        visc_op->maxOrder(max_order);
        Rhs.mult(rhsscale,0,1);
        //
        // Get coefficients for scaled tolerance.
        //
        const MultiFab* alpha = &(visc_op->aCoefficients());
        const MultiFab* beta[BL_SPACEDIM];
        for (int d = 0; d < BL_SPACEDIM; d++)
        {
            beta[d] = &(visc_op->bCoefficients(d));
        }
        //
        // Construct solver and call it.
        //
        const Real S_tol = visc_tol;
        const Real S_tol_abs = get_scaled_abs_tol(sigma, &Rhs, a, b, alpha,
                                                  beta, visc_abs_tol);
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
        delete visc_op;
        //
        // Copy into state variable at new time.
        //
        n_comp  = 1;
        n_ghost = 1;
        Copy(U_new,Soln,0,sigma,n_comp,n_ghost);
        //
        // Modify diffusive fluxes here.
        //
        if (do_reflux)
        {
            FArrayBox flux;

            for (MultiFabIterator U_oldmfi(U_old); U_oldmfi.isValid(); ++U_oldmfi)
            {
                DependentMultiFabIterator U_newmfi(U_oldmfi, U_new);

                int i                  = U_oldmfi.index();
                const Box& grd         = U_oldmfi.validbox();
                const FArrayBox& sold  = U_oldmfi();
                const FArrayBox& snew  = U_newmfi();
                const Real mult        = -visc_coef[sigma];

                for (int d = 0; d < BL_SPACEDIM; d++)
                {
                    flux.resize(::surroundingNodes(grd,d), 1);

                    DependentMultiFabIterator areadmfi(U_oldmfi, area[d]);

                    const FArrayBox& aread = areadmfi();

                    FORT_VISCFLUX_CC(sold.dataPtr(sigma),
                                     snew.dataPtr(sigma),
                                     ARLIM(sold.loVect()), ARLIM(sold.hiVect()),
                                     grd.loVect(), grd.hiVect(),
                                     flux.dataPtr(),
                                     ARLIM(flux.loVect()), ARLIM(flux.hiVect()),
                                     aread.dataPtr(),
                                     ARLIM(aread.loVect()),ARLIM(aread.hiVect()),
                                     &dx[d],&mult,&be_cn_theta,&d);
                    //
                    // NOTE: fluxes expected to be in extensive form.
                    //
                    if (level < finest_level)
                        finer->viscflux_reg->CrseInit(flux,flux.box(),d,0,sigma,1,-dt);
		
                    if (level > 0)
                        viscflux_reg->FineAdd(flux,d,i,0,sigma,1,dt);
                }
            }
            if (level < finest_level)
                finer->viscflux_reg->CrseInitFinish();
        }
    }
}

void
Diffusion::diffuse_tensor_velocity (Real       dt,
                                    Real       be_cn_theta,
                                    MultiFab*  rho_half,
                                    MultiFab*  delta_rhs,
                                    MultiFab** betan, 
                                    MultiFab** betanp1)
{
#if (BL_SPACEDIM==3) 
    cout << "Diffusion::diffuse_tensor_velocity : "
         << "not yet implemented for 3-D\n";
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0"
         << " and rerun\n";
    BoxLib::Abort("Diffusion::diffuse_tensor_velocity");
#elif !defined(USE_TENSOR)
    cout << "Diffusion::diffuse_tensor_velocity : \n";
    cout << "USE_TENSOR must be defined at compile time\n";
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0"
         << " and rerun\n";
    BoxLib::Abort("Diffusion::diffuse_tensor_velocity");
#else
    const int finest_level = parent->finestLevel();
    //
    // At this point, S_old has bndry at time N S_new contains GRAD(SU).
    //
    MultiFab& U_old      = caller->get_old_data(State_Type);
    MultiFab& U_new      = caller->get_new_data(State_Type);
    const Box& domain    = caller->Geom().Domain();
    const int* domlo     = domain.loVect();
    const int* domhi     = domain.hiVect();
    const Real* dx       = caller->Geom().CellSize();
    const Real cur_time  = caller->get_state_data(State_Type).curTime();
    const Real prev_time = caller->get_state_data(State_Type).prevTime();
    //
    // U_new now contains the inviscid update of U.
    // This is part of the RHS for the viscous solve.
    //
    MultiFab Rhs(grids,BL_SPACEDIM,0,Fab_allocate);
    Rhs.setVal(0.0);

    MultiFab** tensorflux_old;
    {
        //
        // Set up Rhs.
        //
        int soln_old_grow = 1;
        MultiFab Soln_old(grids,BL_SPACEDIM,soln_old_grow,Fab_allocate);
        const Real a = 0.0;
        const Real b = -(1.0-be_cn_theta)*dt;

        ViscBndry2D visc_bndry;
        DivVis* tensor_op = getTensorOp(a,b,prev_time,visc_bndry,rho_half,betan);
        tensor_op->maxOrder(tensor_max_order);
        //
        // Copy to single-component multifab.  Use Soln as a temporary here.
        //
        Copy(Soln_old,U_old,Xvel,0,BL_SPACEDIM,soln_old_grow);
        tensor_op->apply(Rhs,Soln_old);

        if (do_reflux && (level<finest_level || level>0))
        {
            allocFluxBoxesLevel(tensorflux_old,0,BL_SPACEDIM);
            tensor_op->compFlux(*(tensorflux_old[0]),*(tensorflux_old[1]),
#if(BL_SPACEDIM==3)
                                *(tensorflux_old[2]),
#endif
                                Soln_old);
            for (int dim = 0; dim < BL_SPACEDIM; dim++)
            {
                tensorflux_old[dim]->mult(-(1.0-be_cn_theta),0);
            }
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
                    delta_rhsmfi().mult(dt,comp,1);
                    delta_rhsmfi().mult(volumemfi(),Rhsmfi.validbox(),0,comp,1);
                    Rhsmfi().plus(delta_rhsmfi(),Rhsmfi.validbox(),comp,comp,1);
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
                const int *lo       = bx.loVect();
                const int *hi       = bx.hiVect();
                const int *slo      = sbx.loVect();
                const int *shi      = sbx.hiVect();
                Real *rhs           = Rhsmfi().dataPtr();
                Real *sdat          = U_oldmfi().dataPtr(Xvel);
                const Real *rcendat = rcen.dataPtr();
                const Real coeff    = (1.0-be_cn_theta)*dt;
                const Real *voli    = volumemfi().dataPtr();
                Box vbox            = ::grow(volumemfi.validbox(),volume.nGrow());
                const int *vlo      = vbox.loVect();
                const int *vhi      = vbox.hiVect();
    
                FArrayBox& betax = betanp10mfi();
                DEF_CLIMITS(betax,betax_dat,betax_lo,betax_hi);

                FArrayBox& betay = betanp11mfi();
                DEF_CLIMITS(betay,betay_dat,betay_lo,betay_hi);

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
    MultiFab Soln(grids,BL_SPACEDIM,soln_grow,Fab_allocate);
    Soln.setVal(0.0);
    //
    // Compute guess of solution.
    //
    if (level == 0)
    {
        Copy(Soln,U_old,Xvel,0,BL_SPACEDIM,soln_grow);
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
    const Real b = be_cn_theta*dt;
       
    ViscBndry2D visc_bndry;
    DivVis* tensor_op = getTensorOp(a,b,cur_time,visc_bndry,rho_half,betanp1);
    tensor_op->maxOrder(tensor_max_order);
    const MultiFab* alpha = &(tensor_op->aCoefficients());
    //
    // Construct solver and call it.
    //
    const Real S_tol = visc_tol;
    const Real S_tol_abs = get_scaled_abs_tol(Xvel, &Rhs, a, b, alpha, betan,
                                              betanp1, visc_abs_tol);
    if (use_tensor_cg_solve)
    {
        int use_mg_pre = 0;
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
    Copy(U_new,Soln,0,Xvel,n_comp,n_ghost);
    //
    // Modify diffusive fluxes here.
    //
    if (do_reflux && (level<finest_level || level>0))
    {
        MultiFab** tensorflux;
        allocFluxBoxesLevel(tensorflux,0,BL_SPACEDIM);
        tensor_op->compFlux(*(tensorflux[0]),*(tensorflux[1]),
#if(BL_SPACEDIM==3)
                            *(tensorflux[2]),
#endif
                            Soln);
        for (int dim = 0; dim < BL_SPACEDIM; dim++)
        {
            tensorflux[dim]->mult(-be_cn_theta,0);
            tensorflux[dim]->plus(*(tensorflux_old[dim]),0,BL_SPACEDIM,0);
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
                assert(grids[tensorflux0mfi.index()] == tensorflux0mfi.validbox());

                int i          = tensorflux0mfi.index();
                const Box& grd = tensorflux0mfi.validbox();
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
                FluxRegister& fr = *finer->viscflux_reg;
                fr.CrseInitFinish();
            }
        }
        removeFluxBoxesLevel(tensorflux);
    }
    delete tensor_op;
#endif /*defined(USE_TENSOR)*/
}

void
Diffusion::diffuse_Vsync (MultiFab*  Vsync,
                          Real       dt,
                          Real       be_cn_theta,
                          MultiFab*  rho_half,
                          int        rho_flag,
                          MultiFab** beta)
{
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "Diffusion::diffuse_Vsync\n";
    }

    int allnull, allthere;
    checkBeta(beta, allthere, allnull);

    int constant_viscosity = allnull;
    int use_dv_constant_mu = use_dv_constant_mu_def;

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        use_dv_constant_mu = use_dv_constant_mu && visc_coef[Xvel+i] >= 0.0;
    }
    if (!use_dv_constant_mu && use_dv_constant_mu_def)
    {
        cout << "Diffusion::diffuse_Vsync : must have velocity visc_coefs "
             << ">= 0.0 if use_dv_constant_mu == 1\n";
        BoxLib::Abort("Diffusion::diffuse_Vsync");
    }

    if (constant_viscosity || use_dv_constant_mu)
    {
        diffuse_Vsync_constant_mu(Vsync,dt,be_cn_theta,rho_half,rho_flag);
    }
    else
    {
        diffuse_tensor_Vsync(Vsync,dt,be_cn_theta,rho_half,rho_flag,beta);
    }
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
Diffusion::diffuse_Vsync_constant_mu (MultiFab* Vsync,
                                      Real      dt,
                                      Real      be_cn_theta,
                                      MultiFab* rho_half,
                                      int       rho_flag)
{
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "Diffusion::diffuse_Vsync\n";
    }
    const Real* dx = caller->Geom().CellSize();
    //
    // At this point in time we can only do decoupled scalar
    // so we loop over components.
    //
    for (int comp = 0; comp < BL_SPACEDIM; comp++)
    {
        MultiFab Soln(grids,1,1,Fab_allocate);
        MultiFab Rhs(grids,1,0,Fab_allocate);

        Soln.setVal(0);
        Rhs.copy(*Vsync,comp,0,1);

        Real r_norm = 0.0;
        for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
        {
            r_norm = Max(r_norm,Rhsmfi().norm(0));
        }
        ParallelDescriptor::ReduceRealMax(r_norm);

        if (ParallelDescriptor::IOProcessor())
        {
            cout << "Original max of Vsync " << r_norm << NL;
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

        Copy(*Vsync,Soln,0,comp,1,1);

        Real s_norm = 0.0;
        for (MultiFabIterator Solnmfi(Soln); Solnmfi.isValid(); ++Solnmfi)
        {
            s_norm = Max(s_norm,Solnmfi().norm(0));
        }
        ParallelDescriptor::ReduceRealMax(s_norm);

        if (ParallelDescriptor::IOProcessor())
        {
            cout << "Final max of Vsync " << s_norm << NL;
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
                int i = Vsyncmfi.index();

                const Box& grd = Vsyncmfi.validbox();
                const int* lo = grd.loVect();
                const int* hi = grd.hiVect();

                FArrayBox& u_sync = Vsyncmfi();

                const int* ulo = u_sync.loVect();
                const int* uhi = u_sync.hiVect();

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
Diffusion::diffuse_tensor_Vsync (MultiFab*  Vsync,
                                 Real       dt,
                                 Real       be_cn_theta,
                                 MultiFab*  rho_half,
                                 int        rho_flag,
                                 MultiFab** beta)
{
#if (BL_SPACEDIM==3) 
    cout << "Diffusion::diffuse_tensor_Vsync : "
         << "not yet implemented for 3-D\n";
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0"
         << " and rerun\n";
    BoxLib::Abort("Exiting.");
#elif !defined(USE_TENSOR)
    cout << "Diffusion::diffuse_tensor_Vsync :  "
         << "USE_TENSOR must be defined at compile time\n";
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0"
         << " and rerun\n";
    BoxLib::Abort("Exiting.");
#else
    const int finest_level = parent->finestLevel();
    const Real* dx         = caller->Geom().CellSize();

    MultiFab Soln(grids,BL_SPACEDIM,1,Fab_allocate);
    MultiFab Rhs(grids,BL_SPACEDIM,0,Fab_allocate);

    Soln.setVal(0);
    Rhs.copy(*Vsync,0,0,BL_SPACEDIM);

    Real r_norm = 0.0;
    for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
        r_norm = Max(r_norm,Rhsmfi().norm(0));
    }
    ParallelDescriptor::ReduceRealMax(r_norm);

    if (ParallelDescriptor::IOProcessor())
    {
        cout << "Original max of Vsync " << r_norm << NL;
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

    DivVis* tensor_op = getTensorOp(a,b,rho_half,beta);
    tensor_op->maxOrder(tensor_max_order);
    //
    // Construct solver and call it.
    //
    const Real S_tol      = visc_tol;
    const MultiFab* alpha = &(tensor_op->aCoefficients());
    MultiFab const * betan[BL_SPACEDIM];
    MultiFab const * betanp1[BL_SPACEDIM];
    for (int d = 0; d < BL_SPACEDIM; d++)
    {
        betan[d] = &tensor_op->bCoefficients(d);
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

    Copy(*Vsync,Soln,0,0,BL_SPACEDIM,1);

    Real s_norm = 0.0;
    for (MultiFabIterator Solnmfi(Soln); Solnmfi.isValid(); ++Solnmfi)
    {
        s_norm = Max(s_norm,Solnmfi().norm(0));
    }
    ParallelDescriptor::ReduceRealMax(s_norm);

    if (ParallelDescriptor::IOProcessor())
    {
        cout << "Final max of Vsync " << s_norm << NL;
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
        for (int dim =0; dim <BL_SPACEDIM; dim++)
        {
            tensorflux[0]->mult(-be_cn_theta*dt*dt,0);
        }
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
                assert(grids[tensorflux0mfi.index()] == tensorflux0mfi.validbox()) ;

                int i          = tensorflux0mfi.index();
                const Box& grd = tensorflux0mfi.validbox();
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
                Real one = 1.0;
                viscflux_reg->FineAdd(xflux,0,i,0,sigma,1,one);
                viscflux_reg->FineAdd(yflux,1,i,0,sigma,1,one);
#if (BL_SPACEDIM == 3)
                viscflux_reg->FineAdd(zflux,2,i,0,sigma,1,one);
#endif
            }
        }
        removeFluxBoxesLevel(tensorflux);
    }
    delete tensor_op;
#endif /*defined(USE_TENSOR)*/
}

void
Diffusion::diffuse_Ssync (MultiFab*  Ssync,
                          int        sigma,
                          Real       dt,
                          Real       be_cn_theta,
                          MultiFab*  rho_half,
                          int        rho_flag,
                          int        do_viscsyncflux, 
                          MultiFab** beta,
                          MultiFab*  alpha)
{
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "Diffusion::diffuse_Ssync for scalar " << sigma << NL;
    }

    int allnull, allthere;
    checkBeta(beta, allthere, allnull);

    const Real* dx = caller->Geom().CellSize();

    MultiFab Soln(grids,1,1,Fab_allocate);
    MultiFab Rhs(grids,1,0,Fab_allocate);

    Soln.setVal(0);
    Rhs.copy(*Ssync,sigma,0,1);

    Real r_norm = 0.0;
    for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
        r_norm = Max(r_norm,Rhsmfi().norm(0));
    }
    ParallelDescriptor::ReduceRealMax(r_norm);

    if (ParallelDescriptor::IOProcessor())
    {
        cout << "Original max of Ssync " << r_norm << NL;
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
    Real b       = be_cn_theta*dt;
    if (allnull)
    {
        b *= visc_coef[BL_SPACEDIM+sigma];
    }
    Real rhsscale = 1.0;
    ABecLaplacian* visc_op = getViscOp(BL_SPACEDIM+sigma,a,b,rho_half,rho_flag,&rhsscale,beta,alpha);
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
    const Real S_tol_abs = get_scaled_abs_tol(BL_SPACEDIM+sigma, &Rhs, a, b,
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

    int visc_op_lev = 0;
    visc_op->applyBC(Soln,visc_op_lev);

    Copy(*Ssync,Soln,0,sigma,1,1);

    Real s_norm = 0.0;
    for (MultiFabIterator Solnmfi(Soln); Solnmfi.isValid(); ++Solnmfi)
    {
        s_norm = Max(s_norm,Solnmfi().norm(0));
    }
    ParallelDescriptor::ReduceRealMax(s_norm);

    if (ParallelDescriptor::IOProcessor())
    {
        cout << "Final max of Ssync " << s_norm << NL;
    }

    delete visc_op;

    FArrayBox xflux, yflux;

    if (level > 0 && do_viscsyncflux == 1)
    {
        for (MultiFabIterator Ssyncmfi(*Ssync); Ssyncmfi.isValid(); ++Ssyncmfi)
        {
            DependentMultiFabIterator area0mfi(Ssyncmfi, area[0]);
            DependentMultiFabIterator area1mfi(Ssyncmfi, area[1]);
#if (BL_SPACEDIM == 3)
            DependentMultiFabIterator area2mfi(Ssyncmfi, area[2]);
#endif
            assert(grids[Ssyncmfi.index()] == Ssyncmfi.validbox());

            const int i       = Ssyncmfi.index();
            const Box& grd    = Ssyncmfi.validbox();
            const int* lo     = grd.loVect();
            const int* hi     = grd.hiVect();
            FArrayBox& s_sync = Ssyncmfi();
            const int* slo    = s_sync.loVect();
            const int* shi    = s_sync.hiVect();

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

            Real mult = -b;

#if (BL_SPACEDIM == 2)
            FORT_VISCSYNCFLUX (s_sync.dataPtr(sigma), ARLIM(slo), ARLIM(shi),
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

            FORT_VISCSYNCFLUX (s_sync.dataPtr(sigma), ARLIM(slo), ARLIM(shi),
                               lo,hi,
                               xflux_dat,ARLIM(xflux_lo),ARLIM(xflux_hi),
                               yflux_dat,ARLIM(yflux_lo),ARLIM(yflux_hi),
                               zflux_dat,ARLIM(zflux_lo),ARLIM(zflux_hi),
                               xarea_dat,ARLIM(xarea_lo),ARLIM(xarea_hi),
                               yarea_dat,ARLIM(yarea_lo),ARLIM(yarea_hi),
                               zarea_dat,ARLIM(zarea_lo),ARLIM(zarea_hi),
                               dx,&mult);
#endif
            if (allthere)
            {
                DependentMultiFabIterator beta0mfi(Ssyncmfi,*beta[0]);
                DependentMultiFabIterator beta1mfi(Ssyncmfi,*beta[1]);
                xflux.mult(beta0mfi());
                yflux.mult(beta1mfi());
#if (BL_SPACEDIM == 3)
                DependentMultiFabIterator beta2mfi(Ssyncmfi,*beta[2]);
                zflux.mult(beta2mfi());
#endif
            }

            Real one = 1.0;
            viscflux_reg->FineAdd(xflux,0,i,0,BL_SPACEDIM+sigma,1,one);
            viscflux_reg->FineAdd(yflux,1,i,0,BL_SPACEDIM+sigma,1,one);
#if (BL_SPACEDIM == 3)
            viscflux_reg->FineAdd(zflux,2,i,0,BL_SPACEDIM+sigma,1,one);
#endif
        }
    }

    if (rho_flag == 2)
    {
        //
        // We just solved for S -- what we want is rho*S.
        //
        MultiFab& S_new = caller->get_new_data(State_Type);
        MultiFab Rho_new(grids,1,1,Fab_allocate);
        Copy(Rho_new,S_new,Density,0,1,1);
        for (MultiFabIterator Ssyncmfi(*Ssync); Ssyncmfi.isValid(); ++Ssyncmfi)
        {
            DependentMultiFabIterator Rho_newmfi(Ssyncmfi, Rho_new);

            Ssyncmfi().mult(Rho_newmfi(),Ssyncmfi.validbox(),0,sigma,1);
        }
    }
    //
    // applyBC has put "incorrect" values in the ghost cells
    // outside external Dirichlet boundaries. Reset these to zero
    // so that conservative interpolation works correctly.
    //
    const BCRec& scalbc = caller->get_desc_lst()[State_Type].getBC(sigma+BL_SPACEDIM);
    Box domain          = ::grow(caller->Geom().Domain(),1);
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        if (scalbc.hi(k) == EXT_DIR)
        {
            IntVect smallend = domain.smallEnd();
            smallend.setVal(k,domain.bigEnd(k));
            Box top_strip(smallend,domain.bigEnd(),IntVect::TheCellVector());
            Ssync->setVal(0.0,top_strip,sigma,1,1);
        }

        if (scalbc.lo(k) == EXT_DIR)
        {
            IntVect bigend = domain.bigEnd();
            bigend.setVal(k,domain.smallEnd(k));
            Box bottom_strip(domain.smallEnd(),bigend,IntVect::TheCellVector());
            Ssync->setVal(0.0,bottom_strip,sigma,1,1);
        }
    }
}

#if defined(USE_TENSOR)
DivVis*
Diffusion::getTensorOp (Real         a,
                        Real         b,
                        Real         time,
#if (BL_SPACEDIM==2)  
                        ViscBndry2D& visc_bndry,
#else
                        ViscBndry3D& visc_bndry,
#endif
                        MultiFab*    rho_half,
                        MultiFab**   beta)
{
#if (BL_SPACEDIM==3) 
    cout << "Diffusion::getTensorOp :  "
         << "not yet implemented for 3-D\n";
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0"
         << " and rerun\n";
    BoxLib::Abort("Diffusion::getTensorOp");
#elif !defined(USE_TENSOR)
    cout << "Diffusion::getTensorOp :  " <<
    cout << "USE_TENSOR must be defined at compile time\n";
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun\n";
    BoxLib::Abort("Diffusion::getTensorOp");
#else
    int allthere;
    checkBeta(beta, allthere);

    const Real* dx = caller->Geom().CellSize();

    getTensorBndryData(visc_bndry,time);

    DivVis* tensor_op = new DivVis(visc_bndry,dx);
    tensor_op->maxOrder(tensor_max_order);

    int isrz   = CoordSys::IsRZ();
    int nghost = 1; // Just like Bill.
    //
    // alpha should be the same size as volume.
    //
    MultiFab alpha(grids,BL_SPACEDIM,nghost,Fab_allocate);
    alpha.setVal(0.0,nghost);

    if (a != 0.0)
    {
        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator volumemfi(alphamfi, volume);
            DependentMultiFabIterator rho_halfmfi(alphamfi, (*rho_half));
            DependentMultiFabIterator beta0mfi(alphamfi, (*beta[0]));
            DependentMultiFabIterator beta1mfi(alphamfi, (*beta[1]));
            assert(alpha.box(alphamfi.index()) == alphamfi.validbox());
            assert(volume.box(alphamfi.index()) == volumemfi.validbox());

            const Box& bx       = alphamfi.validbox();
            Array<Real> rcen(bx.length(0));
            parent->Geom(level).GetCellLoc(rcen, bx, 0);
            const int *lo       = bx.loVect();
            const int *hi       = bx.hiVect();
            Real* alpha_dat     = alphamfi().dataPtr();
            Box abx             = ::grow(bx,alpha.nGrow());
            const int *alo      = abx.loVect();
            const int *ahi      = abx.hiVect();
            const Real *rcendat = rcen.dataPtr();
            const Real *voli    = volumemfi().dataPtr();
            Box vbox            = ::grow(volumemfi.validbox(),volume.nGrow());
            const int *vlo      = vbox.loVect();
            const int *vhi      = vbox.hiVect();

            FArrayBox& Rh = rho_halfmfi();
            DEF_LIMITS(Rh,rho_dat,rlo,rhi);

            FArrayBox& betax = beta0mfi();
            DEF_CLIMITS(betax,betax_dat,betax_lo,betax_hi);

            FArrayBox& betay = beta1mfi();
            DEF_CLIMITS(betay,betay_dat,betay_lo,betay_hi);

            FORT_SET_TENSOR_ALPHA(alpha_dat, ARLIM(alo), ARLIM(ahi),
                                  lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                                  voli, ARLIM(vlo), ARLIM(vhi),
                                  rho_dat,ARLIM(rlo),ARLIM(rhi),
                                  betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                                  betay_dat,ARLIM(betay_lo),ARLIM(betay_hi),&isrz);
        }
    }
    tensor_op->setScalars(a,b);
    tensor_op->aCoefficients(alpha);

    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        MultiFab bcoeffs(area[n].boxArray(),1,nghost,Fab_allocate);
        bcoeffs.setBndry(0);
        Copy(bcoeffs,area[n],0,0,1,nghost);
        for (MultiFabIterator betanmfi(*beta[n]); betanmfi.isValid(); ++betanmfi)
        {
            DependentMultiFabIterator bcoeffsmfi(betanmfi, bcoeffs);
            bcoeffsmfi().mult(dx[n]);
            bcoeffsmfi().mult(betanmfi());
        }
        tensor_op->bCoefficients(bcoeffs,n);
    }

    return tensor_op;
#endif
}

DivVis*
Diffusion::getTensorOp (Real       a,
                        Real       b,
                        MultiFab*  rho_half,
                        MultiFab** beta)
{
#if (BL_SPACEDIM==3) 
    cout << "Diffusion::getTensorOp :  "
         << "not yet implemented for 3-D\n";
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0"
         << " and rerun\n";
    BoxLib::Abort("Diffusion::getTensorOp");
#elif !defined(USE_TENSOR)
    cout << "Diffusion::getTensorOp :  "
         << "USE_TENSOR must be defined at compile time\n";
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0"
         << " and rerun\n";
    BoxLib::Abort("Diffusion::getTensorOp");
#else
    int allthere = beta != 0;
    if (allthere)
    {
        for (int dim = 0; dim < BL_SPACEDIM; dim++)
        {
            allthere = allthere && beta[dim] != 0;
        }
    }
    if (!allthere)
    {
        cout << "Diffusion::getTensorOp() : all betas must allocated\n";
        cout << "  all 0 or all non-0\n";
        BoxLib::Abort("Diffusion::getTensorOp()");
    }

    const Real* dx    = caller->Geom().CellSize();
    const Box& domain = caller->Geom().Domain();

    Array<BCRec> bcarray(2*BL_SPACEDIM);

    for (int idim = 0; idim < BL_SPACEDIM; idim++)
    {
        bcarray[idim] = caller->get_desc_lst()[State_Type].getBC(Xvel+idim);
        bcarray[idim+BL_SPACEDIM] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
                                          D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    }

    IntVect ref_ratio;
    if (level > 0)
    {
        ref_ratio = parent->refRatio(level-1);
    }
    else
    {
        ref_ratio = IntVect::TheUnitVector();
    }
    ViscBndry2D bndry;
    bndry.define(grids,2*BL_SPACEDIM,caller->Geom());
    bndry.setHomogValues(bcarray, ref_ratio[0]);
    DivVis* tensor_op = new DivVis(bndry,dx);
    tensor_op->maxOrder(tensor_max_order);

    int isrz   = CoordSys::IsRZ();
    int nghost = 1; // Just like Bill.
    //
    // alpha should be the same size as volume.
    //
    MultiFab alpha(grids,BL_SPACEDIM,nghost,Fab_allocate);
    alpha.setVal(0.0);

    if (a != 0.0)
    {
        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator volumemfi(alphamfi, volume);
            DependentMultiFabIterator rho_halfmfi(alphamfi, (*rho_half));
            DependentMultiFabIterator beta0mfi(alphamfi, (*beta[0]));
            DependentMultiFabIterator beta1mfi(alphamfi, (*beta[1]));
            assert(alpha.box(alphamfi.index()) == alphamfi.validbox());
            assert(volume.box(alphamfi.index()) == volumemfi.validbox());

            const Box& bx       = alphamfi.validbox();
            Array<Real> rcen(bx.length(0));
            parent->Geom(level).GetCellLoc(rcen, bx, 0);
            const int *lo       = bx.loVect();
            const int *hi       = bx.hiVect();
            Real *alpha_dat     = alphamfi().dataPtr();
            Box abx             = ::grow(bx,alpha.nGrow());
            const int *alo      = abx.loVect();
            const int *ahi      = abx.hiVect();
            const Real *rcendat = rcen.dataPtr();
            const Real *voli    = volumemfi().dataPtr();
            Box vbox            = ::grow(volumemfi.validbox(),volume.nGrow());
            const int *vlo      = vbox.loVect();
            const int *vhi      = vbox.hiVect();

            FArrayBox& Rh = rho_halfmfi();
            DEF_LIMITS(Rh,rho_dat,rlo,rhi);

            FArrayBox& betax = beta0mfi();
            DEF_CLIMITS(betax,betax_dat,betax_lo,betax_hi);

            FArrayBox& betay = beta1mfi();
            DEF_CLIMITS(betay,betay_dat,betay_lo,betay_hi);

            FORT_SET_TENSOR_ALPHA(alpha_dat, ARLIM(alo), ARLIM(ahi),
                                  lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                                  voli, ARLIM(vlo), ARLIM(vhi),
                                  rho_dat,ARLIM(rlo),ARLIM(rhi),
                                  betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                                  betay_dat,ARLIM(betay_lo),ARLIM(betay_hi),&isrz);
        }
    }
    tensor_op->setScalars(a,b);
    tensor_op->aCoefficients(alpha);

    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        MultiFab bcoeffs(area[n].boxArray(),1,nghost,Fab_allocate);
        bcoeffs.setBndry(0);
        Copy(bcoeffs,area[n],0,0,1,nghost);
        for (MultiFabIterator betanmfi(*beta[n]); betanmfi.isValid(); ++betanmfi)
        {
            DependentMultiFabIterator bcoeffsmfi(betanmfi, bcoeffs);
            bcoeffsmfi().mult(dx[n]);
            bcoeffsmfi().mult(betanmfi());
        }
        tensor_op->bCoefficients(bcoeffs,n);
    }

    return tensor_op;
#endif
}
#endif /*defined(USE_TENSOR)*/

ABecLaplacian*
Diffusion::getViscOp (int        comp,
                      Real       a,
                      Real       b,
                      Real       time,
                      ViscBndry& visc_bndry,
                      MultiFab*  rho_half,
                      int        rho_flag, 
                      Real*      rhsscale,
                      MultiFab** beta,
                      MultiFab*  alpha_in)
{
    int allnull, allthere;
    checkBeta(beta, allthere, allnull);

    const Real* dx = caller->Geom().CellSize();

    getBndryData(visc_bndry,comp,1,time,rho_flag);

    ABecLaplacian* visc_op = new ABecLaplacian(visc_bndry,dx);
    visc_op->maxOrder(max_order);

    int usehoop = ((comp==Xvel) && (CoordSys::IsRZ()));
    int useden  = (rho_flag == 1);
    //
    // alpha should be the same size as volume.
    //
    MultiFab alpha(grids,1,GEOM_GROW,Fab_allocate);

    for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
    {
        DependentMultiFabIterator volumemfi(alphamfi, volume);
        DependentMultiFabIterator rho_halfmfi(alphamfi, (*rho_half));
        assert(alpha.box(alphamfi.index()) == alphamfi.validbox());
        assert(volume.box(alphamfi.index()) == volumemfi.validbox());

        const Box& bx       = alphamfi.validbox();
        Array<Real> rcen(bx.length(0));
        parent->Geom(level).GetCellLoc(rcen, bx, 0);
        const int *lo       = bx.loVect();
        const int *hi       = bx.hiVect();
        Real *dat           = alphamfi().dataPtr();
        Box abx             = ::grow(bx,alpha.nGrow());
        const int *alo      = abx.loVect();
        const int *ahi      = abx.hiVect();
        const Real *rcendat = rcen.dataPtr();
        const Real *voli    = volumemfi().dataPtr();
        Box vbox            = ::grow(volumemfi.validbox(),volume.nGrow());
        const int *vlo      = vbox.loVect();
        const int *vhi      = vbox.hiVect();

        FArrayBox& Rh = rho_halfmfi();
        DEF_LIMITS(Rh,rho_dat,rlo,rhi);

        FORT_SETALPHA(dat, ARLIM(alo), ARLIM(ahi),
                      lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                      voli, ARLIM(vlo), ARLIM(vhi),
                      rho_dat,ARLIM(rlo),ARLIM(rhi),&usehoop,&useden);
    }

    //  visc_op->setScalars(a,dx[0]*b);
    if (rho_flag == 2)
    {
        //
        // Using conservative diffing for rho*T.
        //
        MultiFab& S_new = caller->get_new_data(State_Type);
        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator S_newmfi(alphamfi, S_new);
            assert(grids[alphamfi.index()] == alphamfi.validbox());
            alphamfi().mult(S_newmfi(),alphamfi.validbox(),Density,0,1);
        }
    }
    if (alpha_in != 0)
    {
        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator alpha_inmfi(alphamfi, (*alpha_in));
            assert(grids[alphamfi.index()] == alphamfi.validbox());
            alphamfi().mult(alpha_inmfi(),alphamfi.validbox(),0,0,1);
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
                bcoeffsmfi().mult(betanmfi());
                bcoeffsmfi().mult(dx[n]);
            }
            visc_op->bCoefficients(bcoeffs,n);
        }
    }

    return visc_op;
}

ABecLaplacian*
Diffusion::getViscOp (int        comp,
                      Real       a,
                      Real       b,
                      MultiFab*  rho_half,
                      int        rho_flag,
                      Real*      rhsscale,
                      MultiFab** beta,
                      MultiFab*  alpha_in)
{
    int allnull, allthere;
    checkBeta(beta, allthere, allnull);

    const Real* dx    = caller->Geom().CellSize();
    const Box& domain = caller->Geom().Domain();
    const BCRec& bc   = caller->get_desc_lst()[State_Type].getBC(comp);

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
    MultiFab alpha(grids,1,GEOM_GROW,Fab_allocate);

    for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
    {
        DependentMultiFabIterator volumemfi(alphamfi, volume);
        DependentMultiFabIterator rho_halfmfi(alphamfi, (*rho_half));
        assert(alpha.box(alphamfi.index()) == alphamfi.validbox());
        assert(volume.box(alphamfi.index()) == volumemfi.validbox());

        const Box& bx       = alphamfi.validbox();
        Array<Real> rcen(bx.length(0));
        parent->Geom(level).GetCellLoc(rcen, bx, 0);
        const int *lo       = bx.loVect();
        const int *hi       = bx.hiVect();
        Real *dat           = alphamfi().dataPtr();
        Box abx             = ::grow(bx,alpha.nGrow());
        const int *alo      = abx.loVect();
        const int *ahi      = abx.hiVect();
        const Real *rcendat = rcen.dataPtr();
        const Real *voli    = volumemfi().dataPtr();
        Box vbox            = ::grow(volumemfi.validbox(),volume.nGrow());
        const int *vlo      = vbox.loVect();
        const int *vhi      = vbox.hiVect();

        FArrayBox& Rh = rho_halfmfi();
        DEF_LIMITS(Rh,rho_dat,rlo,rhi);

        FORT_SETALPHA(dat, ARLIM(alo), ARLIM(ahi),
                      lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                      voli, ARLIM(vlo), ARLIM(vhi),
                      rho_dat,ARLIM(rlo),ARLIM(rhi),&usehoop,&useden);
    }

    if (rho_flag == 2)
    {
        //
        // Using conservative diffing for rho*S.
        //
        MultiFab& S_new = caller->get_new_data(State_Type);
        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator S_newmfi(alphamfi, S_new);
            assert(grids[alphamfi.index()] == alphamfi.validbox());
            alphamfi().mult(S_newmfi(),alphamfi.validbox(),Density,0,1);
        }
    }
    if (alpha_in != 0)
    {
        for (MultiFabIterator alphamfi(alpha); alphamfi.isValid(); ++alphamfi)
        {
            DependentMultiFabIterator alpha_inmfi(alphamfi, (*alpha_in));
            assert(grids[alphamfi.index()] == alphamfi.validbox());
            alphamfi().mult(alpha_inmfi(),alphamfi.validbox(),0,0,1);
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
                bcoeffsmfi().mult(betanmfi());
            }
            bcoeffs.mult(dx[n],0,1,0);
            visc_op->bCoefficients(bcoeffs,n);
        }
    }

    return visc_op;
}

void
Diffusion::getViscTerms (MultiFab&  visc_terms,
                         int        src_comp,
                         int        comp, 
                         Real       time,
                         int        rho_flag,
                         MultiFab** beta)
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
    const Real* dx = caller->Geom().CellSize();
    MultiFab& S    = caller->get_data(State_Type,time);

    visc_terms.setVal(0.0,comp-src_comp,1,1);
    //
    // FIXME
    // LinOp classes cannot handle multcomponent MultiFabs yet,
    // construct the components one at a time and copy to visc_terms.
    //
    MultiFab visc_tmp(grids,1,1,Fab_allocate);
    MultiFab s_tmp(grids,1,1,Fab_allocate);

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
                    bcoeffsmfi().mult(betanmfi());
                    bcoeffsmfi().mult(dx[n]);
                }
                visc_op.bCoefficients(bcoeffs,n);
            }
        }
        //
        // Copy to single component multifab for operator classes.
        //
        Copy(s_tmp,S,comp,0,1,1);

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
                Box  vbx       = ::grow(bx,visc_tmp.nGrow());
                Box  sbx       = ::grow(s_tmpmfi.validbox(),s_tmp.nGrow());
                Array<Real> rcen(bx.length(0));
                parent->Geom(level).GetCellLoc(rcen, bx, 0);
                const int *lo       = bx.loVect();
                const int *hi       = bx.hiVect();
                const int *vlo      = vbx.loVect();
                const int *vhi      = vbx.hiVect();
                const int *slo      = sbx.loVect();
                const int *shi      = sbx.hiVect();
                Real *vdat          = visc_tmpmfi().dataPtr();
                Real *sdat          = s_tmpmfi().dataPtr();
                const Real *rcendat = rcen.dataPtr();
                const Real mu       = visc_coef[comp];
                FORT_HOOPSRC(ARLIM(lo), ARLIM(hi),
                             vdat, ARLIM(vlo), ARLIM(vhi),
                             sdat, ARLIM(slo), ARLIM(shi),
                             rcendat, &mu);
            }
        }
#endif
        assert(visc_tmp.nGrow() > 0);
        const BoxArray& ba = visc_tmp.boxArray();
        //
        // THIS IS JUST A HACK - DONT KNOW HOW TO FILL VISC TERMS
        // IN GHOST CELLS OUTSIDE FINE GRIDS.
        //
        for (MultiFabIterator visc_tmpmfi(visc_tmp);
             visc_tmpmfi.isValid(); ++visc_tmpmfi)
        {
            assert(ba[visc_tmpmfi.index()] == visc_tmpmfi.validbox());
            const Box& grd  = visc_tmpmfi.validbox();
            const int* lo   = grd.loVect();
            const int* hi   = grd.hiVect();
            FArrayBox& visc = visc_tmpmfi();
            int ncomp = visc.nComp();
            DEF_LIMITS(visc,vdat,vlo,vhi);
            FORT_VISCEXTRAP(vdat,ARLIM(vlo),ARLIM(vhi),lo,hi,&ncomp);
        }
        //
        // Copy into periodic translates of visc_tmp.
        //
        caller->Geom().FillPeriodicBoundary(visc_tmp,false,false);
        //
        // Copy from valid regions of overlapping grids.
        //
        visc_tmp.FillBoundary();

        Copy(visc_terms,visc_tmp,0,comp-src_comp,1,1);
    }
}

void
Diffusion::getTensorViscTerms (MultiFab&  visc_terms, 
                               Real       time,
                               MultiFab** beta)
{
#if (BL_SPACEDIM==3)
    cout << "Diffusion::getTensorViscTerms :  "
         << "not yet implemented for 3-D\n";
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0"
         << " and rerun\n";
    BoxLib::Abort("Diffusion::getTensorViscTerms");
#elif  !defined (USE_TENSOR)
    cout << "Diffusion::getTensorViscTerms :  "
         << "USE_TENSOR must be defined at compile time\n";
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0"
         << " and rerun\n";
    BoxLib::Abort("Diffusion::getTensorViscTerms");
#else
    int allthere;
    checkBeta(beta, allthere);

    int src_comp = Xvel;
    int ncomp    = visc_terms.nComp();
    if (ncomp < BL_SPACEDIM)
    {
        cout << "Diffusion::getTensorViscTerms : visc_terms must have\n";
        cout << "  at least BL_SPACEDIM components\n";
        BoxLib::Abort("Diffusion::getTensorViscTerms");
    }
    int vel_ncomp = BL_SPACEDIM;
    //
    // Before computing the godunov predicitors we may have to
    // precompute the viscous source terms.  To do this we must
    // construct a Laplacian operator, set the coeficients and apply
    // it to the time N data.  First, however, we must precompute the
    // fine N bndry values.  We will do this for each scalar that diffuses.
    //
    const Real* dx = caller->Geom().CellSize();
    MultiFab& S    = caller->get_data(State_Type,time);
    visc_terms.setVal(0.0,src_comp,BL_SPACEDIM,1);
    const int ngrd = grids.length();
    //
    // FIXME
    // LinOp classes cannot handle multcomponent MultiFabs yet,
    // construct the components one at a time and copy to visc_terms.
    //
    MultiFab visc_tmp(grids,BL_SPACEDIM,1,Fab_allocate);
    MultiFab s_tmp(grids,BL_SPACEDIM,1,Fab_allocate);

    if (is_diffusive[src_comp])
    {
        ViscBndry2D visc_bndry;
        getTensorBndryData(visc_bndry,time);
        //
        // Set up operator and apply to compute viscous terms.
        //
        const Real a = 0.0;
        const Real b = -1.0;

        DivVis tensor_op(visc_bndry,dx);
        tensor_op.maxOrder(tensor_max_order);
        tensor_op.setScalars(a,b);

        int nghost = 1; // like bill
        //
        // alpha should be the same size as volume.
        //
        MultiFab alpha(grids,BL_SPACEDIM,nghost,Fab_allocate);
        alpha.setVal(0.0);
        tensor_op.aCoefficients(alpha);

        for (int n = 0; n < BL_SPACEDIM; n++)
        {
            MultiFab bcoeffs(area[n].boxArray(),1,nghost,Fab_allocate);
            bcoeffs.setBndry(0);
            Copy(bcoeffs,area[n],0,0,1,nghost);
            for (MultiFabIterator bcoeffsmfi(bcoeffs); bcoeffsmfi.isValid();
                 ++bcoeffsmfi)
            {
                DependentMultiFabIterator betanmfi(bcoeffsmfi, (*beta[n]));
                bcoeffsmfi().mult(dx[n]);
                bcoeffsmfi().mult(betanmfi());
            }
            tensor_op.bCoefficients(bcoeffs,n);
        }

        Copy(s_tmp,S,Xvel,0,BL_SPACEDIM,1);

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
                const Box& bx       = visc_tmp.box(k);
                Box  vbx            = ::grow(bx,visc_tmp.nGrow());
                Box  sbx            = ::grow(s_tmp.box(k),s_tmp.nGrow());
                Array<Real> rcen(bx.length(0));
                parent->Geom(level).GetCellLoc(rcen, bx, 0);
                const int *lo       = bx.loVect();
                const int *hi       = bx.hiVect();
                const int *vlo      = vbx.loVect();
                const int *vhi      = vbx.hiVect();
                const int *slo      = sbx.loVect();
                const int *shi      = sbx.hiVect();
                Real *vdat          = visc_tmp[k].dataPtr();
                Real *sdat          = s_tmp[k].dataPtr();
                const Real *rcendat = rcen.dataPtr();

                FArrayBox& betax = (*beta[0])[k];
                DEF_CLIMITS(betax,betax_dat,betax_lo,betax_hi);

                FArrayBox& betay = (*beta[1])[k];
                DEF_CLIMITS(betay,betay_dat,betay_lo,betay_hi);

                FORT_TENSOR_HOOPSRC(&fort_xvel_comp,ARLIM(lo), ARLIM(hi),
                                    vdat, ARLIM(vlo), ARLIM(vhi),
                                    sdat, ARLIM(slo), ARLIM(shi),
                                    rcendat, 
                                    betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                                    betay_dat,ARLIM(betay_lo),ARLIM(betay_hi));
            }
        }
#endif
        assert(visc_tmp.nGrow() > 0);
        const BoxArray& ba = visc_tmp.boxArray();
        //
        // THIS IS JUST A HACK - DONT KNOW HOW TO FILL VISC TERMS
        // IN GHOST CELLS OUTSIDE FINE GRIDS.
        //
        for (MultiFabIterator visc_tmpmfi(visc_tmp);
             visc_tmpmfi.isValid(); ++visc_tmpmfi)
        {
            assert(ba[visc_tmpmfi.index()] == visc_tmpmfi.validbox());
            const Box& grd(visc_tmpmfi.validbox());
            const int* lo = grd.loVect();
            const int* hi = grd.hiVect();
            FArrayBox& visc = visc_tmpmfi();
            DEF_LIMITS(visc,vdat,vlo,vhi);
            FORT_VISCEXTRAP(vdat,ARLIM(vlo),ARLIM(vhi),lo,hi,&vel_ncomp);
        }
        //
        // Copy into periodic translates of visc_tmp.
        //
        caller->Geom().FillPeriodicBoundary(visc_tmp,false,false);
        //
        // Copy from valid regions of overlapping grids.
        //
        visc_tmp.FillBoundary();

        Copy(visc_terms,visc_tmp,0,0,BL_SPACEDIM,1);
    }
#endif /*defined (USE_TENSOR)*/
}

void
Diffusion::getBndryData (ViscBndry& bndry,
                         int        src_comp,
                         int        num_comp,
                         Real       time,
                         int        rho_flag)
{
    assert(num_comp == 1);

    const BCRec& bc = caller->get_desc_lst()[State_Type].getBC(src_comp);
    bndry.define(grids,num_comp,caller->Geom());
    //
    // Get state, rho, and have caller set physical bc's.
    //
    MultiFab& S = caller->get_data(State_Type,time);
    assert(S.boxArray() == grids);
    caller->setPhysBoundaryValues(State_Type,src_comp,num_comp,time);
    if (rho_flag == 2)
    { 
        caller->setPhysBoundaryValues(State_Type,Density,1,time);
    }
    //
    // Fill physical boundary values into grow cells of a tmp multifab
    // passed into bndry.  (COI+c-f+periodic handled inside solver)
    //
    // A MultiFab is a huge amount of space in which to pass along
    // the phys bc's. InterpBndryData needs a more efficient interface.
    //
    const int nGrow = 1;
    MultiFab Stmp(grids, num_comp, nGrow, Fab_allocate);
    const Geometry& geom = parent->Geom(level);
    Box grDomainNoPer = ::grow(geom.Domain(), nGrow);
    //
    // Exclude periodic sides and corners.
    //
    if (geom.isAnyPeriodic())
    {
	for (int i = 0; i < BL_SPACEDIM; i++)
	    if (geom.isPeriodic(i))
		grDomainNoPer.grow(i,-nGrow);
    }

    for (MultiFabIterator Stmpmfi(Stmp); Stmpmfi.isValid(); ++Stmpmfi)    
    {
        DependentMultiFabIterator Smfi(Stmpmfi,S);

	BoxList bl = ::boxDiff(Stmpmfi().box() & grDomainNoPer,
                               Stmpmfi().box() & geom.Domain());

	for (BoxListIterator bli(bl); bli; ++bli)
	{
	    Stmpmfi().copy(Smfi(), bli(), src_comp, bli(), 0, num_comp);

	    if (rho_flag == 2)
	    {
                assert(Smfi().min(bli(), Density) > 0.0);
		Stmpmfi().divide(Smfi(), bli(), Density, 0, num_comp);
	    }
	}
    }
    
    if (level == 0)
    {
	bndry.setBndryValues(Stmp,0,0,num_comp,bc);
    }
    else
    {
	BoxArray cgrids = grids;
	cgrids.coarsen(crse_ratio);
	BndryRegister crse_br(cgrids,0,1,1,num_comp);
	crse_br.setVal(BL_BOGUS); // interp for solvers over ALL c-f brs, need safe data
	coarser->FillBoundary(crse_br,src_comp,0,num_comp,time,rho_flag);
	bndry.setBndryValues(crse_br,0,Stmp,0,0,num_comp,crse_ratio,bc);
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

    MultiFab S(caller->boxArray(),num_comp,nGrow,Fab_allocate);
    
    for (FillPatchIterator S_fpi(*caller,S,nGrow,time,State_Type,0,num_comp);
         S_fpi.isValid();
         ++S_fpi)
    {
        S[S_fpi.index()].copy(S_fpi(), 0, 0, num_comp);
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
    bdry.copyFrom(S, nGrow, 0, dest_comp, num_comp);
}

#if defined (USE_TENSOR)
void
Diffusion::getTensorBndryData(
#if (BL_SPACEDIM==2) 
                            ViscBndry2D& bndry, 
#else 
                            ViscBndry3D& bndry, 
#endif
                            Real time)
{
#if (BL_SPACEDIM==3)
    cout << "Diffusion::getTensorBndryData :  "
         << "not yet implemented for 3-D" << NL;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0"
         << " and rerun" << NL;
    BoxLib::Abort("Diffusion::getTensorBndryData()");
#else
    int num_comp = BL_SPACEDIM;
    int src_comp = Xvel;
    //
    // Create the BCRec's interpreted by ViscBndry objects
    //
    Array<BCRec> bcarray(2*BL_SPACEDIM);
    for (int idim = 0; idim < BL_SPACEDIM; idim++)
    {
        bcarray[idim] = caller->get_desc_lst()[State_Type].getBC(src_comp+idim);
        bcarray[idim+BL_SPACEDIM] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
                                          D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    }

    bndry.define(grids,2*num_comp,caller->Geom());

    MultiFab& S = caller->get_data(State_Type,time);
    S.FillBoundary(src_comp,num_comp);
    caller->setPhysBoundaryValues(State_Type,src_comp,num_comp,time);

    if (level == 0)
    {
        bndry.setBndryValues(S,src_comp,0,num_comp,bcarray);
    }
    else
    {
        BoxArray cgrids(grids);
        cgrids.coarsen(crse_ratio);
        BndryRegister crse_br(cgrids,0,1,1,num_comp);
        crse_br.setVal(BL_BOGUS);
        const int rho_flag = 0;
        coarser->FillBoundary(crse_br,src_comp,0,num_comp,time,rho_flag);
        bndry.setBndryValues(crse_br,0,S,src_comp,0,num_comp,crse_ratio[0],bcarray);
    }
}
#endif
#endif /*defined (USE_TENSOR)*/

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
    {
        cout << "Diffusion::checkBetas : all betas must either be"
             << "  all 0 or all non-0\n";
        BoxLib::Abort("Diffusion::checkBetas()");
    }
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
        for (int dim = 0; dim < BL_SPACEDIM; dim++)
        {
            allnull = allnull && beta[dim] == 0;
            allthere = allthere && beta[dim] != 0;
        }
    }
    if (!(allthere || allnull))
    {
        cout << "Diffusion::checkBeta : all betas must either be"
             << "  all 0 or all non-0\n";
        BoxLib::Abort("Diffusion::checkBeta()");
    }
}

void
Diffusion::checkBeta (const MultiFab* const* beta,
                      int&                   allthere) const
{
    allthere = beta != 0;
    if (allthere)
    {
        for (int dim = 0; dim < BL_SPACEDIM; dim++)
        {
            allthere = allthere && beta[dim] != 0;
        }
    }
    if (!allthere)
    {
        cout << "Diffusion::checkBeta : all betas must be"
             << "  all non-0\n";
        BoxLib::Abort("Diffusion::checkBeta()");
    }
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
        fluxbox[dir] = new MultiFab(edge_boxes,nvar,nghost,Fab_allocate);
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
        const int nGrowDU      = 1;
        const Real* dx         = caller->Geom().CellSize();
        NavierStokes& ns_level = *(NavierStokes*) &(parent->getLevel(level));
        MultiFab* divu_fp      = ns_level.getDivCond(nGrowDU,time);

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
Diffusion::compute_divmusi (Real       time,
                            MultiFab** beta,
                            MultiFab&  divmusi)
{
    const int nGrowDU      = 1;
    const Real* dx         = caller->Geom().CellSize();
    NavierStokes& ns_level = *(NavierStokes*) &(parent->getLevel(level));
    MultiFab* divu_fp      = ns_level.getDivCond(nGrowDU,time);

    for (MultiFabIterator divmusimfi(divmusi); divmusimfi.isValid(); ++divmusimfi)
    {
        DependentMultiFabIterator divumfi(divmusimfi,*divu_fp);

        int i           = divmusimfi.index();
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
