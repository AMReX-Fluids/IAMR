//
// $Id: Projection.cpp,v 1.146 2003-02-20 19:18:29 car Exp $
//
#include <winstd.H>

#include <CoordSys.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <MacOpProjDrivers.H>
#include <NavierStokes.H>
#include <Profiler.H>
#include <Projection.H>
#include <PROJECTION_F.H>
#include <NAVIERSTOKES_F.H>
#include <hg_projector.H>
#include "ProjOutFlowBC.H"

#ifndef _NavierStokes_H_
enum StateType {State_Type=0, Press_Type};
#if (BL_SPACEDIM == 2)
enum StateNames  { Xvel=0, Yvel, Density};
#else
enum StateNames  { Xvel=0, Yvel, Zvel, Density};
#endif
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define DEF_BOX_LIMITS(box,boxlo,boxhi)   \
const int* boxlo = (box).loVect();           \
const int* boxhi = (box).hiVect();

#define BogusValue 1.e200
#define MAX_LEV 10

//
// Initialization of static members.
//
int       Projection::verbose            = 0;
int       Projection::P_code             = 0;
Real      Projection::proj_tol           = 1.0e-12;
Real      Projection::sync_tol           = 1.0e-8;
Real      Projection::proj_abs_tol       = 1.0e-16;
Real      Projection::filter_factor      = 0.1;
Real      Projection::divu_minus_s_factor= 0.0;
int       Projection::rho_wgt_vel_proj   = 0;
int       Projection::do_outflow_bcs     = 1;
int       Projection::make_sync_solvable = 0;
int       Projection::proj_0             = 0;
int       Projection::proj_2             = 0;
int       Projection::add_vort_proj      = 0;

#ifndef BL_USE_HGPROJ_SERIAL
#if BL_SPACEDIM == 2
#if BL_PRVERSION == 9
static holy_grail_amr_multigrid::stencil hg_stencil = holy_grail_amr_multigrid::full;
#else
static holy_grail_amr_multigrid::stencil hg_stencil = holy_grail_amr_multigrid::cross;
#endif
#elif BL_SPACEDIM == 3
static holy_grail_amr_multigrid::stencil hg_stencil = holy_grail_amr_multigrid::cross;
#endif
#endif // BL_USE_HGPROJ_SERIAL

static RegType project_bc [] =
{
    interior, inflow, outflow, refWall, refWall, refWall
};

#define LEVEL_PROJ      1001
#define INITIAL_VEL     1002
#define INITIAL_PRESS   1003
#define INITIAL_SYNC    1004

Projection::Projection (Amr*   _parent,
                        BCRec* _phys_bc, 
                        int    _do_sync_proj,
                        int    _finest_level, 
                        int    _radius_grow )
   :
    parent(_parent),
    phys_bc(_phys_bc), 
    do_sync_proj(_do_sync_proj),
    finest_level(_finest_level),
    radius_grow(_radius_grow), 
    LevelData(_finest_level+1), 
    radius(_finest_level+1)
{
    read_params();

    if (verbose && ParallelDescriptor::IOProcessor()) 
        std::cout << "Creating projector\n";

    projector_bndry = 0;

    setUpBcs();

    sync_proj = 0;
}

Projection::~Projection ()
{
    if (verbose && ParallelDescriptor::IOProcessor()) 
        std::cout << "Deleting projector\n";

    delete sync_proj;
    delete projector_bndry;

    sync_proj       = 0;
    projector_bndry = 0;
}

void
Projection::read_params ()
{
    //
    // Read parameters from input file and command line.
    //
    ParmParse pp("proj");

    pp.query("v",verbose);
    pp.query("Pcode",P_code);

    pp.query("proj_tol",proj_tol);
    pp.query("sync_tol",sync_tol);
    pp.query("proj_abs_tol",proj_abs_tol);

    pp.query("make_sync_solvable",make_sync_solvable);

    pp.query("filter_factor",filter_factor);

    pp.query("proj_0",proj_0);
    pp.query("proj_2",proj_2);

    if (proj_0 && proj_2)
	BoxLib::Error("You can only set proj_0 OR proj_2");

    if (!proj_2) 
	BoxLib::Error("With new gravity and outflow stuff, must use proj_2");

    pp.query("divu_minus_s_factor",divu_minus_s_factor);

    pp.query("rho_wgt_vel_proj",rho_wgt_vel_proj);

    pp.query("do_outflow_bcs",do_outflow_bcs);

    pp.query("add_vort_proj",add_vort_proj);

    {
#ifndef BL_USE_HGPROJ_SERIAL
	std::string stencil = "cross";
	if ( pp.query("stencil", stencil) )
	{
	    if ( stencil == "cross" )
	    {
		hg_stencil = holy_grail_amr_multigrid::cross;
	    }
	    else if ( stencil == "terrain" )
	    {
		hg_stencil = holy_grail_amr_multigrid::terrain;
	    }
	    else if ( stencil == "full" )
	    {
		hg_stencil = holy_grail_amr_multigrid::full;
	    }
	    else
	    {
		BoxLib::Error("stencil must be cross, terrain, or full");
	    }
	}
#endif
    }
}

void 
Projection::setUpBcs ()
{
    delete projector_bndry;
    //
    // Set up projector bndry.
    //
    const Geometry& geom  = parent->Geom(0);
    const int*      lo_bc = phys_bc->lo();
    const int*      hi_bc = phys_bc->hi();

    RegType proj_bc[BL_SPACEDIM][2];

    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
        proj_bc[i][0] = project_bc[lo_bc[i]];
        proj_bc[i][1] = project_bc[hi_bc[i]];
        if (geom.isPeriodic(i))
        {
            proj_bc[i][0] = periodic;
            proj_bc[i][1] = periodic;
        }
    }

#ifdef BL_USE_HGPROJ_SERIAL
    projector_bndry = new inviscid_fluid_boundary_class(proj_bc);
#else
    projector_bndry = new inviscid_fluid_boundary(proj_bc);
#endif
}

//
// Install a level of the projection.
//

void
Projection::install_level (int           level,
                           AmrLevel*     level_data,
                           PArray<Real>* _radius)
{
    if (verbose && ParallelDescriptor::IOProcessor()) 
        std::cout << "Installing projector level " << level << '\n';

    finest_level = parent->finestLevel();

    if (level > LevelData.size() - 1) 
    {
        LevelData.resize(finest_level+1);
        radius.resize(finest_level+1);
    }

    LevelData.clear(level);
    LevelData.set(level, level_data);
    radius.clear(level);
    radius.set(level, _radius);

    delete sync_proj;

    sync_proj = 0;
}

//
// Build the aliasLib projection object.
//

void
Projection::bldSyncProject ()
{
    const Box& fdomain = parent->Geom(finest_level).Domain();

    delete sync_proj;

    const Array<IntVect>& ref_ratio = parent->refRatio();
    //
    // Build mesh and ratio arrays for entire hierarchy.
    //
    Array<BoxArray> amesh(finest_level+1);
    Array<IntVect> gen_ratio(finest_level);
    for (int lev = 0; lev <= finest_level; lev++) 
    {
        amesh.set(lev, parent->boxArray(lev));
        if (lev > 0)
            gen_ratio.set(lev-1, ref_ratio[lev-1]);
    }

    if (verbose && ParallelDescriptor::IOProcessor()) 
    {
        std::cout << "bldSyncProject:: amr_mesh = \n";
        amr_multigrid::mesh_write(amesh, gen_ratio, fdomain, std::cout);
    }

#ifdef BL_USE_HGPROJ_SERIAL
#if BL_SPACEDIM == 2
    if (CoordSys::IsRZ())
        amr_multigrid::SetRZ();
#endif
    sync_proj = new holy_grail_amr_projector(amesh, gen_ratio, fdomain,
                                             0, finest_level, finest_level,
                                             *projector_bndry, P_code);
#else
    sync_proj = new holy_grail_amr_projector(amesh, gen_ratio, fdomain,
                                             0, finest_level, finest_level,
                                             *projector_bndry,
					     hg_stencil,
                                             P_code);
#if BL_SPACEDIM == 2
    if (CoordSys::IsRZ())
        sync_proj->setCoordSys(holy_grail_amr_multigrid::rz);
#endif
#endif

#ifdef ATMOSPHERE
    //
    // This is not the usual way of setting parameters.
    //
    sync_proj->line_solve_dim = BL_SPACEDIM - 1;  
#endif
  
    if (make_sync_solvable) 
        sync_proj->make_it_so();
}

//
//  Perform a level projection in the advance function
//  Explanation of arguments to the level projector:
//
//  rho_half  contains rho^{n+1/2}
//  U_old  contains the u^n velocities
//  U_new  starts as u^*, is converted to (u^* - u^n)/dt,
//         becomes (u^{n+1} - u^n)/dt in the solver,
//         and is converted back to u^n+1 at the end
//  P_old  contains p^{n-1/2}
//  P_new  gets cleared, initialized to an intial guess for p^{n+1/2}
//         using coarse grid data if available,
//         becomes pressure update phi in the solver,
//         and then converted into final prssure p^{n+1/2}
//

void
Projection::level_project (int             level,
                           Real            time,
                           Real            dt,
                           Real            cur_pres_time,
                           Real            prev_pres_time,
                           const Geometry& geom, 
                           MultiFab&       U_old,
                           MultiFab&       U_new,
                           MultiFab&       P_old,
                           MultiFab&       P_new,
                           MultiFab*       rho_half, 
                           SyncRegister*   crse_sync_reg, 
                           SyncRegister*   fine_sync_reg,  
                           int             crse_dt_ratio,
                           int**           bc,
                           int             iteration,
                           int             have_divu, 
                           int             Divu_Type) 
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::level_project()");

    if (ParallelDescriptor::IOProcessor() && verbose)
	std::cout << "... level projector at level " << level << '\n';

    if (sync_proj == 0)
        bldSyncProject();
    //
    // old time velocity has bndry values already
    // must gen valid bndry data for new time velocity.
    // must fill bndry cells in pressure with computable values
    // even though they are not used in calculation.
    //
    U_old.setBndry(BogusValue,Xvel,BL_SPACEDIM);
    U_new.setBndry(BogusValue,Xvel,BL_SPACEDIM);
    P_old.setBndry(BogusValue);
    P_new.setBndry(BogusValue);
    LevelData[level].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,0);
    LevelData[level].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,1);

    const Real*     dx      = geom.CellSize();
    const BoxArray& grids   = LevelData[level].boxArray();
    const BoxArray& P_grids = P_old.boxArray();

    NavierStokes* ns = dynamic_cast<NavierStokes*>(&parent->getLevel(level));
    BL_ASSERT(!(ns==0));

    //
    //  NOTE: IT IS IMPORTANT TO DO THE BOUNDARY CONDITIONS BEFORE
    //    MAKING UNEW HOLD U_t OR U/dt, BECAUSE UNEW IS USED IN
    //    CONSTRUCTING THE OUTFLOW BC'S.
    //
    // Set boundary values for P_new, to increment, if applicable
    //
    // Note: we don't need to worry here about using FillCoarsePatch because
    //       it will automatically use the "new dpdt" to interpolate,
    //       since once we get here at level > 0, we've already defined
    //       a new pressure at level-1.
    if (level != 0)
    {
	LevelData[level].FillCoarsePatch(P_new,0,cur_pres_time,Press_Type,0,1);
        if (!proj_2) 
            P_new.minus(P_old,0,1,0); // Care about nodes on box boundary
    }
    const int nGrow = (level == 0  ?  0  :  -1);
    for (MFIter P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi)
    {
        const int i = P_newmfi.index();

        P_new[i].setVal(0.0,BoxLib::grow(P_new.box(i),nGrow),0,1);
        //
        // TODO -- also zero fine-fine nodes ???
        //
    }

    //
    // Compute Ustar/dt as input to projector for proj_0,
    //         Ustar/dt + Gp                  for proj_2,
    //         (Ustar-Un)/dt for not proj_0 or proj_2 (ie the original).
    //
    // Compute DU/dt for proj_0,
    //         DU/dt for proj_2,
    //         (DU-DU_old)/dt for not proj_0 or proj_2 (ie the original).
    //
    MultiFab *divusource = 0, *divuold = 0;

    if (have_divu)
    {
        divusource = ns->getDivCond(1,time+dt);
        if (!proj_0 && !proj_2)
            divuold = ns->getDivCond(1,time);
    }

    const Real dt_inv = 1./dt;
    if (proj_0 || proj_2)
    {
        U_new.mult(dt_inv,0,BL_SPACEDIM,1);
        if (have_divu)
            divusource->mult(dt_inv,0,1,divusource->nGrow());
    }
    else
    {
        for (MFIter U_newmfi(U_new); U_newmfi.isValid(); ++U_newmfi) 
        {
            const int i = U_newmfi.index();

            ConvertUnew(U_new[i],U_old[i],dt,U_new.box(i));
        }
        if (have_divu)
        {
            divusource->minus(*divuold,0,1,divusource->nGrow());
            divusource->mult(dt_inv,0,1,divusource->nGrow());

            if (divu_minus_s_factor>0.0 && divu_minus_s_factor<=1.0)
            {
                BoxLib::Error("Check this code....not recently tested");
                //
                // Compute relaxation terms to account for approximate projection
                // add divu_old*divu...factor/dt to divusource.
                //
                const Real uoldfactor = divu_minus_s_factor*dt/parent->dtLevel(0);
                UpdateArg1(*divusource, uoldfactor/dt, *divuold, 1, grids, 1);
                //
                // add U_old*divu...factor/dt to U_new
                //
                UpdateArg1(U_new, uoldfactor/dt, U_old, BL_SPACEDIM, grids, 1);
            }
        }
    }
    delete divuold;

    if (proj_2)
    {
        MultiFab Gp(grids,BL_SPACEDIM,1);

        ns->getGradP(Gp, prev_pres_time);

        const int n_ghost = 1;

        rho_half->invert(1.0,n_ghost);

        for (FillPatchIterator P_fpi(LevelData[level],P_old,1,prev_pres_time,Press_Type,0,1);
             P_fpi.isValid();
             ++P_fpi)
        {
            const int idx = P_fpi.index();

            for (int i = 0; i < BL_SPACEDIM; i++)
                Gp[idx].mult((*rho_half)[idx],0,i,1);

            U_new[idx].plus(Gp[idx],0,0,BL_SPACEDIM);
        }

        rho_half->invert(1.0,n_ghost);
    }

    //
    // Outflow uses appropriately constructed "U_new" and "divusource"
    //   so make sure this call comes after those are set,
    //   but before fields are scaled by r or rho is set to 1/rho.
    //
    Real gravity = ns->getGravity();
    if (OutFlowBC::HasOutFlowBC(phys_bc) && (have_divu || gravity > 0.0) 
                                         && do_outflow_bcs) 
    {
        MultiFab* phi[MAX_LEV] = {0};
        phi[level] = &LevelData[level].get_new_data(Press_Type);

        MultiFab* Vel_ML[MAX_LEV] = {0};
        Vel_ML[level] = &U_new;

        MultiFab* Divu_ML[MAX_LEV] = {0};
        Divu_ML[level] = divusource;

        MultiFab* Rho_ML[MAX_LEV] = {0};
        Rho_ML[level] = rho_half;

        set_outflow_bcs(LEVEL_PROJ,phi,Vel_ML,Divu_ML,Rho_ML,level,level,have_divu);
    }

    //
    // Scale the projection variables.
    //
    rho_half->setBndry(BogusValue);
    scaleVar(rho_half, 1, &U_new, grids, level);

    //
    // Application specific first guess.
    //
    ProjFirstGuess(U_new, P_new, level, grids);
    //
    // Enforce periodicity of U_new and rho_half (i.e. coefficient of G phi)
    // *after* everything has been done to them.
    //
    EnforcePeriodicity(U_new,     BL_SPACEDIM, grids, geom);
    EnforcePeriodicity(*rho_half, 1,           grids, geom);
    //
    // Add the contribution from the un-projected V to syncregisters.
    //
    int is_rz = (CoordSys::IsRZ() ? 1 : 0);
    //
    // Setup projection (note that u_real is a temporary copy).
    //
    // Build amr_real data structures that projection wants.
    //
    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> p_real(level+1), s_real(level+1);
    int n;
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        u_real[n].resize(level+1,PArrayManage);
        u_real[n].set(level, new MultiFab(grids, 1, 1));
        for (MFIter u_realmfi(u_real[n][level]); u_realmfi.isValid(); ++u_realmfi)
        {
            const int i = u_realmfi.index();

            u_real[n][level][i].copy(U_new[i], n, 0);
        }
    }
    p_real.set(level, &P_new);
    s_real.set(level, rho_half);
    //
    // Project
    //
    MultiFab* sync_resid_crse = 0;
    MultiFab* sync_resid_fine = 0;

    if (level < finest_level) 
        sync_resid_crse = new MultiFab(P_grids,1,1);

    if ( level > 0 &&
         ( ((proj_0 || proj_2) && iteration == crse_dt_ratio) ||
           (!proj_0 && !proj_2)) )
    {
        const int ngrow = parent->MaxRefRatio(level-1) - 1;
        sync_resid_fine = new MultiFab(P_grids,1,ngrow);
    }

    if (!have_divu) 
    {
        sync_proj->project(u_real, p_real, null_amr_real, s_real, 
                           sync_resid_crse, sync_resid_fine, geom, 
                           (Real*)dx, proj_tol, level, level, proj_abs_tol);
    }
    else 
    {
        bool use_u = true;
        if (is_rz == 1)
            radMult(level,*divusource,0);
        const int nghost = 0;
        divusource->mult(-1.0,0,1,nghost); // FIXME: this doesn't touch the ghost cells?

        PArray<MultiFab> rhs_real(level+1);
        rhs_real.set(level, divusource);
        sync_proj->manual_project(u_real,p_real,rhs_real,null_amr_real,s_real, 
                                  sync_resid_crse, sync_resid_fine, geom, 
                                  use_u, (Real*)dx,
                                  proj_tol, level, level, proj_abs_tol);
    }

    delete divusource;
    //
    // Note: this must occur *after* the projection has been done
    //       (but before the modified velocity has been copied back)
    //       because the SyncRegister routines assume the projection
    //       has been set up.
    //
    if (do_sync_proj)
    {
       if (level < finest_level)
       {
          //
          // Init sync registers between level and level+1.
          //
          const Real mult = 1.0;
          crse_sync_reg->CrseInit(sync_resid_crse,geom,mult);
       }
       if ( level > 0 &&
            ( ((proj_0 || proj_2) && iteration == crse_dt_ratio) ||
              (!proj_0 && !proj_2)) )
       {
          //
          // Increment sync registers between level and level-1.
          //
          const Real      invrat    = 1.0/(double)crse_dt_ratio;
          const Geometry& crse_geom = parent->Geom(level-1);
          fine_sync_reg->FineAdd(sync_resid_fine,geom,crse_geom,phys_bc,invrat);
       }
    }
    delete sync_resid_crse;
    delete sync_resid_fine;
    //
    // Copy u_real.
    //
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        for (MFIter u_realmfi(u_real[n][level]); u_realmfi.isValid(); ++u_realmfi)
        {
            const int i = u_realmfi.index();

            U_new[i].copy(u_real[n][level][i], 0, n);
        }
    }
    //
    // Reset state + pressure data.
    //
    // Unscale level projection variables.
    //
    rescaleVar(rho_half, 1, &U_new, grids, level);
    //
    // Put U_new back to "normal"; subtract U_old*divu...factor/dt from U_new
    //
    if ( (!proj_0 && !proj_2) && 
        divu_minus_s_factor>0.0 && divu_minus_s_factor<=1.0 && have_divu) 
    {
        const Real uoldfactor = -divu_minus_s_factor*dt/parent->dtLevel(0);
        UpdateArg1(U_new, uoldfactor/dt, U_old, BL_SPACEDIM, grids, 1);
    }
    //
    // Convert U back to a velocity, and phi into p^n+1/2.
    //
    if (proj_0 || proj_2) 
    {
        //
        // un = dt*un
        //
        U_new.mult(dt,0,BL_SPACEDIM,1);
    }
    else
    {
        //
        // un = uo+dt*un
        //
        UnConvertUnew(U_old, dt, U_new, grids);
    }

    if (!proj_2) 
        AddPhi(P_new, P_old, grids);             // pn = pn + po

    if (proj_0)
    {
        const Real dt_inv = 1./dt;
        U_old.mult(-dt_inv,0,BL_SPACEDIM,1);
        filterP(level,geom,P_old,P_new,U_old,rho_half,bc,time,dt,have_divu);
        U_old.mult(-dt,0,BL_SPACEDIM,1);
    }
}

void
Projection::filterP (int             level,
                     const Geometry& geom, 
                     MultiFab&       P_old,
                     MultiFab&       P_new,
                     MultiFab&       U_old,
                     MultiFab*       rho_half, 
                     int**           bc,
                     Real            time,
                     Real            dt,
                     int             have_divu)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... filterP at level " << level << std::endl;

    int             rzflag   = CoordSys::IsRZ();
    const Real*     dx       = geom.CellSize();
    Box             ndomain  = BoxLib::surroundingNodes(geom.Domain());
    const int*      ndlo     = ndomain.loVect();
    const int*      ndhi     = ndomain.hiVect();
    const BoxArray& grids    = LevelData[level].boxArray();
    const BoxArray& P_grids  = P_old.boxArray();
    const int*      phys_lo  = phys_bc->lo();
    const int*      phys_hi  = phys_bc->hi();
    MultiFab*       rhs_nd   = new MultiFab(P_grids,1,0);
    MultiFab*       temp_phi = new MultiFab(P_grids,1,1);
    MultiFab*       temp_rho = new MultiFab(grids,1,1);
    MultiFab*       temp_vel = new MultiFab(grids,BL_SPACEDIM,1);
    MultiFab*       sync_resid_crse = 0;
    MultiFab*       sync_resid_fine = 0;

    BL_ASSERT(grids.size() == P_grids.size());
    
    temp_phi->setVal(0);
    temp_rho->setVal(0);
    rhs_nd->setVal(0);
    //
    // Scale the projection variables.
    //
    scaleVar(rho_half, 0, &U_old, grids, level);
    //
    // Copy from valid regions only.
    //
    temp_rho->copy(*rho_half,0,0,1);
    temp_rho->FillBoundary();
    
    temp_phi->copy(P_old,0,0,1);
    temp_phi->FillBoundary();
    //
    // Copy from valid regions + bcs to get inflow values.
    //
    int n_ghost = 1;
    MultiFab::Copy(*temp_vel,U_old,0,0,BL_SPACEDIM,n_ghost);

    EnforcePeriodicity(*temp_vel, BL_SPACEDIM, grids, geom);
    EnforcePeriodicity(*temp_rho, 1,           grids, geom);
    EnforcePeriodicity(*temp_phi, 1,           P_grids, geom);

    Real mult = -1.;

    for (MFIter mfi(*temp_rho); mfi.isValid(); ++mfi)
    {
        const int  k    = mfi.index();
        FArrayBox& sfab = (*temp_rho)[k];
        const int* s_lo = sfab.loVect();
        const int* s_hi = sfab.hiVect();
        FArrayBox& pfab = (*temp_phi)[k];
        const int* p_lo = pfab.loVect();
        const int* p_hi = pfab.hiVect();
        const int* r_lo = (*rhs_nd)[k].loVect();
        const int* r_hi = (*rhs_nd)[k].hiVect();
        const int* n_lo = P_grids[k].loVect();
        const int* n_hi = P_grids[k].hiVect();

        FORT_FILTRHS(pfab.dataPtr(),ARLIM(p_lo),ARLIM(p_hi),
                     sfab.dataPtr(),ARLIM(s_lo),ARLIM(s_hi),
                     (*rhs_nd)[k].dataPtr(),ARLIM(r_lo),ARLIM(r_hi),
                     n_lo,n_hi,dx,&mult,&rzflag);

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (!geom.isPeriodic(dir))
            {
                Box ndomlo(ndomain), ndomhi(ndomain);
                ndomlo.setRange(dir,ndlo[dir],1);
                ndomhi.setRange(dir,ndhi[dir],1);
                //
                // Any RHS point on the physical bndry must be multiplied by
                // two (only for ref-wall, inflow and symmetry) and set to 
                // zero at outflow.
                //
                Box blo = P_grids[k] & ndomlo;
                Box bhi = P_grids[k] & ndomhi;

                if (blo.ok())
                {
                    if (phys_lo[dir] == Outflow) 
                        (*rhs_nd)[k].setVal(0,blo,0,1);
                    else
                        (*rhs_nd)[k].mult(2,blo,0,1);
                }
                if (bhi.ok())
                {
                    if (phys_hi[dir] == Outflow) 
                        (*rhs_nd)[k].setVal(0,bhi,0,1);
                    else
                        (*rhs_nd)[k].mult(2,bhi,0,1);
                }
            } 
        }
    }

    if (have_divu)
    {
        const int nghost  = 0;

        MultiFab* divuold = new MultiFab(P_grids,1,nghost);

        put_divu_in_node_rhs(*divuold,level,nghost,time,rzflag);

        for (MFIter mfi(*divuold); mfi.isValid(); ++mfi)
        {
            (*divuold)[mfi].mult(1.0/dt,0,1);
        }

        rhs_nd->plus(*divuold,0,1,nghost);

        delete divuold;
    }

    temp_phi->setVal(0);
    //
    // Setup projection .
    //
    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> p_real(level+1), s_real(level+1);
    PArray<MultiFab> rhs_real(level+1);
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        u_real[n].resize(level+1,PArrayManage);
        u_real[n].set(level, new MultiFab(grids, 1, 1));
        for (MFIter u_realmfi(u_real[n][level]); u_realmfi.isValid(); ++u_realmfi)
        {
            const int i = u_realmfi.index();

            u_real[n][level][i].copy((*temp_vel)[i], n, 0);
        }
    }

    p_real.set(level, temp_phi);
    s_real.set(level, temp_rho);
    rhs_real.set(level, rhs_nd);

    delete temp_vel;
    delete temp_rho;
    delete rhs_nd;
    //
    // Project ...
    //
    const bool use_u   = true;

    if (level < finest_level) 
        sync_resid_crse = new MultiFab(P_grids,1,1);

    if (level > 0) 
    {
        const int ngrow = parent->MaxRefRatio(level-1) - 1;
        sync_resid_fine = new MultiFab(P_grids,1,ngrow);
    }

    sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                              sync_resid_crse, sync_resid_fine, geom, 
                              use_u, (Real*)dx,
                              filter_factor, level, level, proj_abs_tol);
    //
    // Reset state + pressure data
    //
    AddPhi(P_new, *temp_phi, grids);
    //
    // Unscale the projection variables.
    //
    rescaleVar(rho_half, 0, &U_old, grids, level);

    delete temp_phi;
    delete sync_resid_crse;
    delete sync_resid_fine;
}

//
// SYNC_PROJECT
//

#define MAXLEV 10

void
Projection::syncProject (int             c_lev,
                         MultiFab&       pres,
                         MultiFab&       vel,
                         MultiFab*       rho_half,
                         MultiFab*       Vsync,
                         MultiFab&       phi,
                         SyncRegister*   rhs_sync_reg,
                         SyncRegister*   crsr_sync_reg,
                         const BoxArray& sync_boxes,
                         const Geometry& geom,
                         const Real*     dx,
                         Real            dt_crse,
                         int             crse_iteration,
                         int             crse_dt_ratio)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::syncProject()");

    int rz_flag = (CoordSys::IsRZ() ? 1 : 0);

    if (verbose && ParallelDescriptor::IOProcessor()) 
    {
        std::cout << "SyncProject: level = "
                  << c_lev
                  << " correction to level "
                  << finest_level << std::endl;
    }
    //
    // Manipulate state + pressure data.
    //
    if (sync_proj == 0)
        bldSyncProject();
    //
    // Gather data.
    //
    const BoxArray& grids   = LevelData[c_lev].boxArray();
    const BoxArray& P_grids = pres.boxArray();
    MultiFab  rhs(P_grids,1,1);
    MultiFab& sig = *rho_half;

    rhs_sync_reg->InitRHS(rhs,geom,phys_bc);

    phi.setVal(0);

    sig.setBndry(BogusValue);
    //
    // Scale sync projection variables.
    //
    scaleVar(&sig,1,Vsync,grids,c_lev);
    //
    // If periodic, copy into periodic translates of Vsync.
    //
    EnforcePeriodicity(*Vsync, BL_SPACEDIM, grids, geom);
    //
    // Build aliaslib structures.
    //
    int n;
    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> p_real(c_lev+1), s_real(c_lev+1), rhs_real(c_lev+1);
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        u_real[n].resize(c_lev+1,PArrayManage);
        u_real[n].set(c_lev, new MultiFab(grids, 1, 1));
        for (MFIter u_realmfi(u_real[n][c_lev]); u_realmfi.isValid(); ++u_realmfi)
        {
            const int i = u_realmfi.index();

            u_real[n][c_lev][i].copy((*Vsync)[i], n, 0);
        }
    }
    p_real.set(c_lev, &phi);
    s_real.set(c_lev, &sig);
    rhs_real.set(c_lev, &rhs);
    //
    //  PROJECT
    //  if use_u = 0, then solves DGphi = RHS
    //  if use_u = 1, then solves DGphi = RHS + DV
    //  both return phi and (V-Gphi) as V
    //
    const bool use_u          = true;
    MultiFab* sync_resid_crse = 0;
    MultiFab* sync_resid_fine = 0;

    if (c_lev > 0)
    {
        const int ngrow = parent->MaxRefRatio(c_lev-1) - 1;
        sync_resid_fine = new MultiFab(P_grids,1,ngrow);
    }

    sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                              sync_resid_crse, sync_resid_fine, geom, 
                              use_u, (Real*)dx,
                              sync_tol, c_lev, c_lev, proj_abs_tol);
    //
    // Copy u_real.
    //
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        for (MFIter u_realmfi(u_real[n][c_lev]); u_realmfi.isValid(); ++u_realmfi)
        {
            const int i = u_realmfi.index();

            (*Vsync)[i].copy(u_real[n][c_lev][i], 0, n);
        }
    }
    //
    // If this sync project is not at level 0 then we need to account for
    // the changes made here in the level c_lev velocity in the sync registers
    // going into the level (c_lev-1) sync project.  Note that this must be
    // done before rho_half is scaled back.
    //
    if (c_lev > 0 &&
        (!proj_2 || crse_iteration == crse_dt_ratio) )
    {
        const Real invrat         = 1.0/(double)crse_dt_ratio;
        const Geometry& crsr_geom = parent->Geom(c_lev-1);
        crsr_sync_reg->CompAdd(sync_resid_fine,geom,crsr_geom,phys_bc,sync_boxes,invrat);
    }
    delete sync_resid_crse;
    delete sync_resid_fine;
    //
    // Reset state + pressure data ...
    //
    // Unscale the sync projection variables for rz.
    //
    rescaleVar(&sig,1,Vsync,grids,c_lev);
    //
    // Add projected Vsync to new velocity at this level & add phi to pressure.
    //
    AddPhi(pres, phi, grids);
    UpdateArg1(vel, dt_crse, *Vsync, BL_SPACEDIM, grids, 1);
}

//
//  MULTI-LEVEL SYNC_PROJECT
//

void
Projection::MLsyncProject (int             c_lev,
                           MultiFab&       pres_crse,
                           MultiFab&       vel_crse,
                           MultiFab&       cc_rhs_crse,
                           MultiFab&       pres_fine,
                           MultiFab&       vel_fine,
                           MultiFab&       cc_rhs_fine,
                           MultiFab&       rho_crse,
                           MultiFab&       rho_fine,
                           MultiFab*       Vsync,
                           MultiFab&       V_corr,
                           MultiFab&       phi_fine,
                           SyncRegister*   rhs_sync_reg,
                           SyncRegister*   crsr_sync_reg,
                           const Real*     dx,
                           Real            dt_crse, 
                           IntVect&        ratio,
                           int             crse_iteration,
                           int             crse_dt_ratio,
                           const Geometry& fine_geom,
                           const Geometry& crse_geom,
                           bool		   pressure_time_is_interval,
                           bool first_crse_step_after_initial_iters,
                           Real             cur_crse_pres_time,
                           Real            prev_crse_pres_time,
                           Real             cur_fine_pres_time,
                           Real            prev_fine_pres_time)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::MLsyncProject()");

    if (verbose && ParallelDescriptor::IOProcessor()) 
        std::cout << "SyncProject: levels = " << c_lev << ", " << c_lev+1 << '\n';
    
    int rz_flag = (CoordSys::IsRZ() ? 1 : 0);
    if (sync_proj == 0)
        bldSyncProject();
    //
    // Set up memory.
    //
    MultiFab *phi[MAXLEV];

    const BoxArray& grids      = LevelData[c_lev].boxArray();
    const BoxArray& fine_grids = LevelData[c_lev+1].boxArray();
    const BoxArray& Pgrids_crse = pres_crse.boxArray();

    phi[c_lev] = new MultiFab(Pgrids_crse,1,1);
    phi[c_lev]->setVal(0);
    
    MultiFab* crse_rhs = new MultiFab(Pgrids_crse,1,1);
    
    const BoxArray& Pgrids_fine = pres_fine.boxArray();
    phi[c_lev+1] = new MultiFab(Pgrids_fine,1,1);
    phi[c_lev+1]->setVal(0);

    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> p_real(c_lev+2,PArrayManage);
    PArray<MultiFab> s_real(c_lev+2), crse_rhs_real(c_lev+1);
    PArray<MultiFab> rhs_real(c_lev+2);
    //
    // Set up crse RHS
    //
    rhs_sync_reg->InitRHS(*crse_rhs,crse_geom,phys_bc);

    Box P_finedomain(BoxLib::surroundingNodes(crse_geom.Domain()));
    P_finedomain.refine(ratio);
    if (Pgrids_fine[0] == P_finedomain)
        crse_rhs->setVal(0);
    //
    // Do necessary scaling
    //
    scaleVar(&rho_crse, 0, Vsync,   grids,      c_lev  );
    scaleVar(&rho_fine, 0, &V_corr, fine_grids, c_lev+1);
    //
    // Set up alias lib.
    //
    int n;
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        u_real[n].resize(c_lev+2,PArrayManage);
        u_real[n].set(c_lev,   new MultiFab(grids,      1, 1));
        u_real[n].set(c_lev+1, new MultiFab(fine_grids, 1, 1));

        for (MFIter u_realmfi(u_real[n][c_lev]); u_realmfi.isValid(); ++u_realmfi)
        {
            const int i = u_realmfi.index();

            u_real[n][c_lev][i].copy((*Vsync)[i], n, 0);
        }
        for (MFIter u_realfinemfi(u_real[n][c_lev+1]); u_realfinemfi.isValid();
             ++u_realfinemfi)
        {
            const int i = u_realfinemfi.index();

            u_real[n][c_lev+1][i].copy(V_corr[i], n, 0);
        }

        restrict_level(u_real[n][c_lev], u_real[n][c_lev+1], ratio);
    }

    s_real.set(c_lev,   &rho_crse);
    s_real.set(c_lev+1, &rho_fine);

    if (CoordSys::IsRZ()) {
       radMult(c_lev  ,cc_rhs_crse,0);
       radMult(c_lev+1,cc_rhs_fine,0);
    }

    rhs_real.set(c_lev  , &cc_rhs_crse);
    rhs_real.set(c_lev+1, &cc_rhs_fine);

    restrict_level(s_real[c_lev], s_real[c_lev+1], ratio);

    p_real.set(c_lev,   phi[c_lev]);
    p_real.set(c_lev+1, phi[c_lev+1]);
    //
    // Note that crse_rhs_real is only built on the coarsest level.
    //
    crse_rhs_real.set(c_lev, crse_rhs);
    //
    // The Multilevel Projection
    // if use_u = false, then solves DGphi = RHS
    // if use_u = true , then solves DGphi = RHS + DV
    // both return phi and (V-Gphi) as V
    //
    const bool  use_u           = true;
    const Real* dx_fine         = parent->Geom(c_lev+1).CellSize();
    MultiFab*   sync_resid_crse = 0;
    MultiFab*   sync_resid_fine = 0;

    if (c_lev > 0)
    {
        int ngrow = parent->MaxRefRatio(c_lev-1) - 1;
        sync_resid_fine = new MultiFab(Pgrids_crse,1,ngrow);
    }
    sync_proj->manual_project(u_real, p_real, rhs_real,
                              crse_rhs_real, s_real, 
                              sync_resid_crse, sync_resid_fine, crse_geom, 
                              use_u, (Real*)dx_fine,
                              sync_tol, c_lev, c_lev+1, proj_abs_tol);
    delete crse_rhs;
    //
    // Copy u_real.
    //
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        for (MFIter u_realmfi(u_real[n][c_lev]); u_realmfi.isValid(); ++u_realmfi)
        {
            const int i = u_realmfi.index();

            (*Vsync)[i].copy(u_real[n][c_lev][i], 0, n);
        }
        for (MFIter u_realfinemfi(u_real[n][c_lev+1]); u_realfinemfi.isValid();
             ++u_realfinemfi)
        {
            const int i = u_realfinemfi.index();

            V_corr[i].copy(u_real[n][c_lev+1][i], 0, n);
        }
    }
    //
    // If this sync project is not at levels 0-1 then we need to account for
    // the changes made here in the level c_lev velocity in the sync registers
    // going into the level (c_lev-1) sync project.  Note that this must be
    // done before rho_half is scaled back.
    //
    if (c_lev > 0 &&
        (!proj_2 || crse_iteration == crse_dt_ratio) )
    {
        const Real invrat         = 1.0/(double)crse_dt_ratio;
        const Geometry& crsr_geom = parent->Geom(c_lev-1);
        BoxArray sync_boxes       = pres_fine.boxArray();
        sync_boxes.coarsen(ratio);
        crsr_sync_reg->CompAdd(sync_resid_fine,crse_geom,crsr_geom,
                               phys_bc,sync_boxes,invrat);
    }
    delete sync_resid_fine;
    //
    // Do necessary un-scaling.
    //
    rescaleVar(&rho_crse, 0, Vsync,   grids,      c_lev  );
    rescaleVar(&rho_fine, 0, &V_corr, fine_grids, c_lev+1);

    for (MFIter phimfi(*phi[c_lev+1]); phimfi.isValid(); ++phimfi) 
    {
        phi_fine[phimfi].copy((*phi[c_lev+1])[phimfi],0,0,1);
    }
    //
    // Add phi to pressure.
    //
    AddPhi(pres_crse, *phi[c_lev],   grids);

    if (pressure_time_is_interval) 
    {
        //
        // Only update the most recent pressure.
        //
        AddPhi(pres_fine, *phi[c_lev+1], fine_grids);
    }
    else 
    {
        MultiFab& pres_fine_old = LevelData[c_lev+1].get_old_data(Press_Type);
 
        if (first_crse_step_after_initial_iters)
        {
            Real time_since_zero =  cur_crse_pres_time - prev_crse_pres_time;
            Real dt_to_prev_time = prev_fine_pres_time - prev_crse_pres_time;
            Real dt_to_cur_time  =  cur_fine_pres_time - prev_crse_pres_time;

            Real cur_mult_factor = dt_to_cur_time / time_since_zero;
            (*phi[c_lev+1]).mult(cur_mult_factor);
            AddPhi(pres_fine, *phi[c_lev+1], fine_grids);

            Real prev_mult_factor = dt_to_prev_time / dt_to_cur_time;
            (*phi[c_lev+1]).mult(prev_mult_factor);
            AddPhi(pres_fine_old, *phi[c_lev+1], fine_grids);
        }
        else 
        {
            AddPhi(pres_fine    , *phi[c_lev+1], fine_grids);
            AddPhi(pres_fine_old, *phi[c_lev+1], fine_grids);
        }
    }
    //
    // Add projected vel to new velocity.
    //
    UpdateArg1(vel_crse, dt_crse, *Vsync, BL_SPACEDIM, grids,      1);
    UpdateArg1(vel_fine, dt_crse, V_corr, BL_SPACEDIM, fine_grids, 1);
}

//
// The initial velocity projection in post_init.
// this function ensures that the velocities are nondivergent
//

void
Projection::initialVelocityProject (int  c_lev,
                                    Real cur_divu_time, 
                                    int  have_divu)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::initialVelocityProject()");

    int lev;
    int f_lev = finest_level;
    if (verbose && ParallelDescriptor::IOProcessor()) 
    {
        std::cout << "initialVelocityProject: levels = " << c_lev
                  << "  " << f_lev << '\n';
        if (rho_wgt_vel_proj) 
            std::cout << "RHO WEIGHTED INITIAL VELOCITY PROJECTION\n";
        else 
            std::cout << "CONSTANT DENSITY INITIAL VELOCITY PROJECTION\n";
    }

    if (sync_proj == 0)
        bldSyncProject();

    MultiFab* vel[MAX_LEV] = {0};
    MultiFab* phi[MAX_LEV] = {0};
    MultiFab* sig[MAX_LEV] = {0};

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        LevelData[lev].get_old_data(Press_Type).setVal(0);
    }

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev] = &LevelData[lev].get_new_data(State_Type);
        phi[lev] = &LevelData[lev].get_old_data(Press_Type);

        const int       nghost = 1;
        const BoxArray& grids  = LevelData[lev].boxArray();
        sig[lev]               = new MultiFab(grids,1,nghost);

        if (rho_wgt_vel_proj) 
        {
            LevelData[lev].get_new_data(State_Type).setBndry(BogusValue,Density,1);

            parent->getLevel(lev).setPhysBoundaryValues(State_Type,Density,1);

            MultiFab::Copy(*sig[lev],
                           LevelData[lev].get_new_data(State_Type),
                           Density,
                           0,
                           1,
                           nghost);
        }
        else 
        {
            sig[lev]->setVal(1,nghost);
        }
    }

    MultiFab* rhs_cc[MAX_LEV] = {0};
    const int nghost = 1; 

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev]->setBndry(BogusValue,Xvel,BL_SPACEDIM);
        //
        // Set the physical boundary values.
        //
        AmrLevel&       amr_level = parent->getLevel(lev);
        const BoxArray& grids     = amr_level.boxArray();

        amr_level.setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM);

        if (have_divu) 
        {
            int Divu_Type, Divu;
            if (!LevelData[lev].isStateVariable("divu", Divu_Type, Divu)) 
                BoxLib::Error("Projection::initialVelocityProject(): Divu not found");
            //
            // Make sure ghost cells are properly filled.
            //
            MultiFab& divu_new = amr_level.get_new_data(Divu_Type);
            divu_new.FillBoundary();
            amr_level.setPhysBoundaryValues(Divu_Type,0,1,cur_divu_time);

            const BoxArray& grids     = amr_level.boxArray();
            rhs_cc[lev]  = new MultiFab(grids,1,nghost);
            MultiFab* rhslev = rhs_cc[lev];
            put_divu_in_cc_rhs(*rhslev,lev,grids,cur_divu_time);
        }
    }

    if (OutFlowBC::HasOutFlowBC(phys_bc) && do_outflow_bcs)
    {
       if (have_divu)
         set_outflow_bcs(INITIAL_VEL,phi,vel,rhs_cc,sig,
                        c_lev,f_lev,have_divu);

       if (!have_divu) {
          MultiFab* divu_dummy[MAX_LEV] = {0};
          set_outflow_bcs(INITIAL_VEL,phi,vel,divu_dummy,sig,
                        c_lev,f_lev,have_divu);
       }
    }

     //
     // Scale the projection variables.
     //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
       AmrLevel&       amr_level = parent->getLevel(lev);
       const BoxArray& grids     = amr_level.boxArray();
       scaleVar(sig[lev],1,vel[lev],grids,lev);
    }

    //
    // Setup alias lib.
    //
    const Array<BoxArray>& full_mesh = sync_proj->mesh();
    PArray<MultiFab> u_real[BL_SPACEDIM], p_real(f_lev+1);
    PArray<MultiFab> s_real(f_lev+1,PArrayManage);

    int n;
    for (n = 0; n < BL_SPACEDIM; n++) 
        u_real[n].resize(f_lev+1,PArrayManage);

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            u_real[n].set(lev, new MultiFab(full_mesh[lev], 1, 1));

            for (MFIter u_realmfi(u_real[n][lev]); u_realmfi.isValid(); ++u_realmfi)
            {
                const int i = u_realmfi.index();

                u_real[n][lev][i].copy((*vel[lev])[i], Xvel+n, 0);
            }
        }
        p_real.set(lev, phi[lev]);
        s_real.set(lev, sig[lev]);
    }
    //
    // Project
    //
    MultiFab*   sync_resid_crse = 0;
    MultiFab*   sync_resid_fine = 0;
    const Real* dx_lev          = parent->Geom(f_lev).CellSize();

    if (!have_divu)
    {
        sync_proj->project(u_real, p_real, null_amr_real, s_real,
                           sync_resid_crse, sync_resid_fine, parent->Geom(c_lev), 
                           (Real*)dx_lev, proj_tol, c_lev, f_lev,proj_abs_tol);
    } 
    else 
    {
        const bool use_u = true;

        PArray<MultiFab> rhs_real(f_lev+1,PArrayManage);
        for (lev = c_lev; lev <= f_lev; lev++) 
        {
            MultiFab* rhslev = rhs_cc[lev];
            if (CoordSys::IsRZ()) radMult(lev,*rhslev,0); 
            rhs_cc[lev]->mult(-1.0,0,1,nghost);
            rhs_real.set(lev, rhs_cc[lev]);
        }

        sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real, 
                                  sync_resid_crse, sync_resid_fine, 
                                  parent->Geom(c_lev), 
                                  use_u, (Real*)dx_lev,
                                  proj_tol, c_lev, f_lev, proj_abs_tol);
    }

    //
    // Copy u_real.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            for (MFIter u_realmfi(u_real[n][lev]); u_realmfi.isValid();
                 ++u_realmfi)
            {
                const int i = u_realmfi.index();

                (*vel[lev])[i].copy(u_real[n][lev][i], 0, Xvel+n);
            }
        }
    }
    //
    // Unscale initial projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        const BoxArray& grids = parent->getLevel(lev).boxArray();
        rescaleVar(sig[lev],1,vel[lev],grids,lev);
    }

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        LevelData[lev].get_old_data(Press_Type).setVal(0);
        LevelData[lev].get_new_data(Press_Type).setVal(0);
    }
}

void
Projection::initialPressureProject (int  c_lev)
                                    
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::initialPressureProject()");

    int lev;
    int f_lev = finest_level;
    if (verbose && ParallelDescriptor::IOProcessor()) 
        std::cout << "initialPressureProject: levels = " << c_lev
                  << "  " << f_lev << '\n';

    MultiFab* vel[MAX_LEV] = {0};
    MultiFab* phi[MAX_LEV] = {0};
    MultiFab* sig[MAX_LEV] = {0};

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev] = &LevelData[lev].get_new_data(State_Type);
        phi[lev] = &LevelData[lev].get_old_data(Press_Type);

        const int       nghost = 1;
        const BoxArray& grids  = LevelData[lev].boxArray();
        sig[lev]               = new MultiFab(grids,1,nghost);

        LevelData[lev].get_new_data(State_Type).setBndry(BogusValue,Density,1);

        parent->getLevel(lev).setPhysBoundaryValues(State_Type,Density,1);

        MultiFab::Copy(*sig[lev],
                       LevelData[lev].get_new_data(State_Type),
                       Density,
                       0,
                       1,
                       nghost);
    }

    //
    // Set up outflow bcs.
    //
    NavierStokes* ns = dynamic_cast<NavierStokes*>(&LevelData[c_lev]);
    Real gravity = ns->getGravity();

    if (OutFlowBC::HasOutFlowBC(phys_bc) && do_outflow_bcs)
    {
        int have_divu_dummy = 0;
        MultiFab* Divu_ML[MAX_LEV] = {0};

        set_outflow_bcs(INITIAL_PRESS,phi,vel,Divu_ML,sig,
                        c_lev,f_lev,have_divu_dummy);
    }

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        AmrLevel&       amr_level = parent->getLevel(lev);
        const BoxArray& grids     = amr_level.boxArray();

        //
        // Scale the projection variables.
        //
        scaleVar(sig[lev],1,vel[lev],grids,lev);
    }

    //
    // Setup alias lib.
    //
    const Array<BoxArray>& full_mesh = sync_proj->mesh();
    PArray<MultiFab> p_real(f_lev+1);
    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> s_real(f_lev+1,PArrayManage);

    int n;
    for (n = 0; n < BL_SPACEDIM; n++) 
        u_real[n].resize(f_lev+1,PArrayManage);

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            u_real[n].set(lev, new MultiFab(full_mesh[lev], 1, 1));

            for (MFIter u_realmfi(u_real[n][lev]); u_realmfi.isValid(); ++u_realmfi)
            {
                const int i = u_realmfi.index();
                if (n == (BL_SPACEDIM-1)) {
                  u_real[n][lev][i].setVal(-gravity);
                } else { 
                  u_real[n][lev][i].setVal(0.);
                }
            }
        }
        p_real.set(lev, phi[lev]);
        s_real.set(lev, sig[lev]);
    }
    //
    // Project
    //
    const Real* dx_lev          = parent->Geom(f_lev).CellSize();
    MultiFab*   sync_resid_crse = 0;
    MultiFab*   sync_resid_fine = 0;

    sync_proj->project(u_real, p_real, null_amr_real, s_real,
                       sync_resid_crse, sync_resid_fine, parent->Geom(c_lev), 
                       (Real*)dx_lev, proj_tol, c_lev, f_lev,proj_abs_tol);

    //
    // Unscale initial projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        const BoxArray& grids = parent->getLevel(lev).boxArray();
        rescaleVar(sig[lev],1,vel[lev],grids,lev);
    }

    //
    // Copy "old" pressure just computed into "new" pressure as well.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
        MultiFab::Copy(LevelData[lev].get_new_data(Press_Type),
                       LevelData[lev].get_old_data(Press_Type),
                       0, 0, 1, 0);
}

//
// The velocity projection in post_init, which computes the initial
// pressure used in the timestepping.
//

void
Projection::initialSyncProject (int       c_lev,
                                MultiFab* sig[],
                                Real      dt, 
                                Real      strt_time,
                                Real      dt_init,
                                int       have_divu)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::initialSyncProject()");

    int lev;
    int f_lev = finest_level;
    if (verbose && ParallelDescriptor::IOProcessor()) 
        std::cout << "initialSyncProject: levels = " << c_lev << "  " << f_lev << '\n';
    //
    // Manipulate state + pressure data.
    //
    if (sync_proj == 0)
        bldSyncProject();
    //
    // Gather data.
    //
    MultiFab* vel[MAX_LEV] = {0};
    MultiFab* phi[MAX_LEV] = {0};
    MultiFab* rhs[MAX_LEV] = {0};

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev] = &LevelData[lev].get_new_data(State_Type);
        phi[lev] = &LevelData[lev].get_old_data(Press_Type);
    }
  
    const Real dt_inv = 1./dt;

    if (have_divu) 
    {
        //
        // Set up rhs for manual project.
        //
        for (lev = c_lev; lev <= f_lev; lev++) 
        {
            AmrLevel& amr_level = parent->getLevel(lev);

            int Divu_Type, Divu;
            if (!LevelData[c_lev].isStateVariable("divu", Divu_Type, Divu)) 
                BoxLib::Error("Projection::initialSyncProject(): Divu not found");
            //
            // Make sure ghost cells are properly filled.
            //
            MultiFab& divu_new = amr_level.get_new_data(Divu_Type);
            divu_new.FillBoundary();
            MultiFab& divu_old = amr_level.get_old_data(Divu_Type);
            divu_old.FillBoundary();
            amr_level.setPhysBoundaryValues(Divu_Type,0,1,strt_time);
            amr_level.setPhysBoundaryValues(Divu_Type,0,1,strt_time+dt);

            const int nghost = 1;
            rhs[lev] = new MultiFab(amr_level.boxArray(),1,nghost);
            MultiFab* rhslev = rhs[lev];
            rhslev->setVal(0);

            NavierStokes* ns = dynamic_cast<NavierStokes*>(&parent->getLevel(lev));

            BL_ASSERT(!(ns == 0));

            MultiFab* divu = ns->getDivCond(nghost,strt_time);
            MultiFab* dsdt = ns->getDivCond(nghost,strt_time+dt);

            for (MFIter mfi(*rhslev); mfi.isValid(); ++mfi)
            {

                (*dsdt)[mfi].minus((*divu)[mfi]);
                (*dsdt)[mfi].mult(dt_inv);
                (*rhslev)[mfi].copy((*dsdt)[mfi]);
            }

            delete divu;
            delete dsdt;
        }
    }

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        MultiFab& P_old = LevelData[lev].get_old_data(Press_Type);
        P_old.setVal(0);
    }
    //
    // Set velocity bndry values to bogus values.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev]->setBndry(BogusValue,Xvel,BL_SPACEDIM);
        MultiFab &u_o = LevelData[lev].get_old_data(State_Type);
        u_o.setBndry(BogusValue,Xvel,BL_SPACEDIM);
        sig[lev]->setBndry(BogusValue);
    }
    //
    // Convert velocities to accelerations (we always do this for the
    //  projections in these initial iterations).
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        LevelData[lev].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,1);
        LevelData[lev].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,0);
        MultiFab& u_o = LevelData[lev].get_old_data(State_Type);
        ConvertUnew(*vel[lev], u_o, dt, LevelData[lev].boxArray());
    }

    if (OutFlowBC::HasOutFlowBC(phys_bc) && have_divu && do_outflow_bcs) 
        set_outflow_bcs(INITIAL_SYNC,phi,vel,rhs,sig,
                        c_lev,f_lev,have_divu);

    //
    // Scale initial sync projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        AmrLevel& amr_level = parent->getLevel(lev);
        scaleVar(sig[lev],1,vel[lev],amr_level.boxArray(),lev);

        if (have_divu && CoordSys::IsRZ()) 
          radMult(lev,*(rhs[lev]),0);    
    }

    //
    // Build the aliaslib data structures
    //
    const Array<BoxArray>& full_mesh = sync_proj->mesh();
    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> p_real(f_lev+1), s_real(f_lev+1);
    int n;
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        u_real[n].resize(f_lev+1,PArrayManage);
    }
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            u_real[n].set(lev, new MultiFab(full_mesh[lev], 1, 1));
            for (MFIter u_realmfi(u_real[n][lev]); u_realmfi.isValid(); ++u_realmfi)
            {
                const int i = u_realmfi.index();

                u_real[n][lev][i].copy((*vel[lev])[i], Xvel+n, 0);
            }
        }
        p_real.set(lev, phi[lev]);
        s_real.set(lev, sig[lev]);
    }

    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        for (lev = f_lev; lev >= c_lev+1; lev--) 
        {
            restrict_level(u_real[n][lev-1],u_real[n][lev],
                           parent->refRatio(lev-1));
        }
    }

    MultiFab*   sync_resid_crse = 0;
    MultiFab*   sync_resid_fine = 0;
    const Real* dx_lev          = parent->Geom(f_lev).CellSize();
    //
    // Project.
    //
    if (!have_divu) 
    {
        //
        // Zero divu only or debugging.
        //
        sync_proj->project(u_real, p_real, null_amr_real, s_real,
                           sync_resid_crse, sync_resid_fine, 
                           parent->Geom(c_lev), 
                           (Real*)dx_lev, proj_tol, c_lev, f_lev,proj_abs_tol);
    } 
    else 
    {
        //
        // General divu.
        //
        bool use_u = true;
        //
        // This PArray is managed so it'll remove rhs on exiting scope.
        //
        PArray<MultiFab> rhs_real(f_lev+1,PArrayManage);
        for (lev = c_lev; lev <= f_lev; lev++) 
        {
            rhs[lev]->mult(-1.0,0,1);
            rhs_real.set(lev, rhs[lev]);
        }
        sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                                  sync_resid_crse, sync_resid_fine, 
                                  parent->Geom(c_lev), 
                                  use_u, (Real*)dx_lev,
                                  proj_tol, c_lev, f_lev, proj_abs_tol);
    }
    //
    // Copy u_real.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            for (MFIter u_realmfi(u_real[n][lev]); u_realmfi.isValid(); ++u_realmfi)
            {
                const int i = u_realmfi.index();

                (*vel[lev])[i].copy(u_real[n][lev][i], 0, Xvel+n);
            }
        }
    }
    //
    // Unscale initial sync projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
        rescaleVar(sig[lev],1,vel[lev],parent->getLevel(lev).boxArray(),lev);
    //
    // Add correction at coarse and fine levels.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
        incrPress(lev, 1.0);
}

//
// Put S in the rhs of the projector -- node based version.
//

void
Projection::put_divu_in_node_rhs (MultiFab&       rhs,
                                  int             level,
                                  const int&      nghost,
                                  Real            time,
                                  int             user_rz)
{
    BL_ASSERT(user_rz >= -1 && user_rz <= 1);

    rhs.setVal(0);

    const Geometry& geom   = parent->Geom(level);
    const int       bcxlo  = phys_bc->lo(0);
    const int       bcxhi  = phys_bc->hi(0);
    int             isrz   = user_rz == -1 ? CoordSys::IsRZ() : user_rz;
    const Real*     dx     = geom.CellSize();
    Real            hx     = dx[0];
    const Box&      domain = geom.Domain();
    const int*      domlo  = domain.loVect();
    const int*      domhi  = domain.hiVect();

    NavierStokes* ns = dynamic_cast<NavierStokes*>(&parent->getLevel(level));

    BL_ASSERT(!(ns == 0));

    MultiFab* divu = ns->getDivCond(1,time);
 
    int imax = geom.Domain().bigEnd()[0]+1;

    for (MFIter rhsmfi(rhs); rhsmfi.isValid(); ++rhsmfi)
    {
        DEF_CLIMITS((*divu)[rhsmfi],divudat,divulo,divuhi);
        DEF_LIMITS(rhs[rhsmfi],rhsdat,rhslo,rhshi);
        Array<Real> rcen((*divu)[rhsmfi].box().length(0),1.0);
        if (isrz == 1) 
            geom.GetCellLoc(rcen,(*divu)[rhsmfi].box(),0);

        FORT_HGC2N(&nghost,ARLIM(divulo),ARLIM(divuhi),divudat,
                   rcen.dataPtr(), ARLIM(rhslo),ARLIM(rhshi),rhsdat,
                   domlo,domhi,&hx,&isrz,&imax);
    }

    delete divu;
}

//
// Put S in the rhs of the projector--cell based version.
//

void
Projection::put_divu_in_cc_rhs (MultiFab&       rhs,
                                int             level,
                                const BoxArray& grids,
                                Real            time)
{
    rhs.setVal(0);

    NavierStokes* ns = dynamic_cast<NavierStokes*>(&parent->getLevel(level));

    BL_ASSERT(!(ns == 0));

    MultiFab* divu = ns->getDivCond(1,time);

    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
        rhs[mfi].copy((*divu)[mfi]);
    }

    delete divu;
}

void
Projection::EnforcePeriodicity (MultiFab&       psi,
                                int             nvar,
                                const BoxArray& /*grids*/,
                                const Geometry& geom)
{
    BL_ASSERT(nvar <= psi.nComp());

    geom.FillPeriodicBoundary(psi,0,nvar);
}

//
// Convert U from an Accl-like quantity to a velocity: Unew = Uold + alpha*Unew
//

void
Projection::UnConvertUnew (MultiFab&       Uold,
                           Real            alpha,
                           MultiFab&       Unew, 
                           const BoxArray& grids)
{
    for (MFIter Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        BL_ASSERT(grids[Uoldmfi.index()] == Uoldmfi.validbox());

        UnConvertUnew(Uold[Uoldmfi],alpha,Unew[Uoldmfi],Uoldmfi.validbox());
    }
}

//
// Convert U from an Accleration like quantity to a velocity
// Unew = Uold + alpha*Unew.
//

void
Projection::UnConvertUnew (FArrayBox& Uold,
                           Real       alpha,
                           FArrayBox& Unew,
                           const Box& grd)
{
    BL_ASSERT(Unew.nComp() >= BL_SPACEDIM);
    BL_ASSERT(Uold.nComp() >= BL_SPACEDIM);
    BL_ASSERT(Unew.contains(grd) == true);
    BL_ASSERT(Uold.contains(grd) == true);
    
    const int*  lo    = grd.loVect();
    const int*  hi    = grd.hiVect();
    const int*  uo_lo = Uold.loVect(); 
    const int*  uo_hi = Uold.hiVect(); 
    const Real* uold  = Uold.dataPtr(0);
    const int*  un_lo = Unew.loVect(); 
    const int*  un_hi = Unew.hiVect(); 
    const Real* unew  = Unew.dataPtr(0);
    
    FORT_ACCEL_TO_VEL(lo, hi,
                      uold, ARLIM(uo_lo), ARLIM(uo_hi),
                      &alpha,
                      unew, ARLIM(un_lo), ARLIM(un_hi));
}

//
// Convert U to an Accleration like quantity: Unew = (Unew - Uold)/alpha
//

void
Projection::ConvertUnew (MultiFab&       Unew,
                         MultiFab&       Uold,
                         Real            alpha,
                         const BoxArray& grids)
{
    for (MFIter Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        BL_ASSERT(grids[Uoldmfi.index()] == Uoldmfi.validbox());

        ConvertUnew(Unew[Uoldmfi],Uold[Uoldmfi],alpha,Uoldmfi.validbox());
    }
}

//
// Convert U to an Accleration like quantity: Unew = (Unew - Uold)/alpha
//

void
Projection::ConvertUnew( FArrayBox &Unew, FArrayBox &Uold, Real alpha,
                              const Box &grd )
{
    BL_ASSERT(Unew.nComp() >= BL_SPACEDIM);
    BL_ASSERT(Uold.nComp() >= BL_SPACEDIM);
    BL_ASSERT(Unew.contains(grd) == true);
    BL_ASSERT(Uold.contains(grd) == true);
    
    const int*  lo    = grd.loVect();
    const int*  hi    = grd.hiVect();
    const int*  uo_lo = Uold.loVect(); 
    const int*  uo_hi = Uold.hiVect(); 
    const Real* uold  = Uold.dataPtr(0);
    const int*  un_lo = Unew.loVect(); 
    const int*  un_hi = Unew.hiVect(); 
    const Real* unew  = Unew.dataPtr(0);
                    
    FORT_VEL_TO_ACCEL(lo, hi, 
                      unew, ARLIM(un_lo), ARLIM(un_hi),
                      uold, ARLIM(uo_lo), ARLIM(uo_hi), &alpha );
}

//
// Update a quantity U using the formula: Unew = Unew + alpha*Uold
//

void
Projection::UpdateArg1 (MultiFab&       Unew,
                        Real            alpha,
                        MultiFab&       Uold,
                        int             nvar,
                        const BoxArray& grids,
                        int             ngrow)
{
    for (MFIter Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        BL_ASSERT(grids[Uoldmfi.index()] == Uoldmfi.validbox());

        UpdateArg1(Unew[Uoldmfi],alpha,Uold[Uoldmfi],nvar,Uoldmfi.validbox(),ngrow);
    }
}

//
// Update a quantity U using the formula
// currently only the velocity, but will do the pressure as well.
// Unew = Unew + alpha*Uold
//

void
Projection::UpdateArg1 (FArrayBox& Unew,
                        Real       alpha,
                        FArrayBox& Uold,
                        int        nvar,
                        const Box& grd,
                        int        ngrow)
{
    BL_ASSERT(nvar <= Uold.nComp());
    BL_ASSERT(nvar <= Unew.nComp());

    Box        b  = BoxLib::grow(grd,ngrow);
    const Box& bb = Unew.box();

    if (bb.ixType() == IndexType::TheNodeType())
        b.surroundingNodes();

    BL_ASSERT(Uold.contains(b) == true);
    BL_ASSERT(Unew.contains(b) == true);

    const int*  lo    = b.loVect();
    const int*  hi    = b.hiVect();
    const int*  uo_lo = Uold.loVect(); 
    const int*  uo_hi = Uold.hiVect(); 
    const Real* uold  = Uold.dataPtr(0);
    const int*  un_lo = Unew.loVect(); 
    const int*  un_hi = Unew.hiVect(); 
    const Real* unew  = Unew.dataPtr(0);
                    
    FORT_PROJ_UPDATE(lo,hi,&nvar,&ngrow,
                     unew, ARLIM(un_lo), ARLIM(un_hi),
                     &alpha,
                     uold, ARLIM(uo_lo), ARLIM(uo_hi) );
}

//
// Add phi to P.
//

void
Projection::AddPhi (MultiFab&        p,
                    MultiFab&       phi,
                    const BoxArray& grids)
{
    for (MFIter pmfi(p); pmfi.isValid(); ++pmfi) 
    {
        p[pmfi].plus(phi[pmfi]);
    }
}

//
// Convert phi into p^n+1/2.
//

void
Projection::incrPress (int  level,
                       Real dt)
{
    MultiFab& P_old = LevelData[level].get_old_data(Press_Type);
    MultiFab& P_new = LevelData[level].get_new_data(Press_Type);

    const BoxArray& grids = LevelData[level].boxArray();

    for (MFIter P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi)
    {
        const int i = P_newmfi.index();

        UpdateArg1(P_new[P_newmfi],1.0/dt,P_old[P_newmfi],1,grids[i],1);

        P_old[P_newmfi].setVal(BogusValue);
    }
}

//
// This function scales variables at the start of a projection.
//

void
Projection::scaleVar (MultiFab*       sig,
                      int             sig_nghosts,
                      MultiFab*       vel,
                      const BoxArray& grids,
                      int             level)
{
    if (sig != 0)
        BL_ASSERT(sig->nComp() == 1);
    if (vel != 0)
        BL_ASSERT(vel->nComp() >= BL_SPACEDIM);
    //
    // Convert sigma from rho to 1/rho.
    // nghosts info needed to avoid divide by zero.
    //
    if (sig != 0)
        sig->invert(1.0,sig_nghosts);
    //
    // Scale by radius for RZ.
    //
    if (CoordSys::IsRZ()) 
    {
        if (sig != 0)
            radMult(level,*sig,0);
        if (vel != 0)
            for (int n = 0; n < BL_SPACEDIM; n++) 
                radMult(level,*vel,n);
    }
    //
    // Scale level projection variables for a particular projection.
    //
    proj_scale_var(sig,vel,grids,level);
}

//
// This function rescales variables at the end of a projection.
//

void
Projection::rescaleVar (MultiFab*       sig,
                        int             sig_nghosts,
                        MultiFab*       vel,
                        const BoxArray& grids,
                        int             level)
{
    if (sig != 0)
        BL_ASSERT(sig->nComp() == 1);
    if (vel != 0)
        BL_ASSERT(vel->nComp() >= BL_SPACEDIM);
    //
    // Divide by radius to rescale for RZ coordinates.
    //
    if (CoordSys::IsRZ()) 
    {
        if (sig != 0)
            radDiv(level,*sig,0);
        if (vel != 0)
            for (int n = 0; n < BL_SPACEDIM; n++)
                radDiv(level,*vel,n);
    }
    //
    // Convert sigma from 1/rho to rho
    // NOTE: this must come after division by r to be correct,
    // nghosts info needed to avoid divide by zero.
    //
    if (sig != 0)
        sig->invert(1.0,sig_nghosts);
    //
    // Unscale level projection variables for a particular projection.
    //
    proj_unscale_var(sig,vel,grids,level);
}

//
// Multiply by a radius for r-z coordinates.
//

void
Projection::radMult (int       level,
                     MultiFab& mf,
                     int       comp)
{
    BL_ASSERT(radius_grow >= mf.nGrow());
    BL_ASSERT(comp >= 0 && comp < mf.nComp());

    int ngrow = mf.nGrow();

    int nr = radius_grow;

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

    for (MFIter mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
    {
        BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

        const int* lo = mfmfi.validbox().loVect();
        const int* hi = mfmfi.validbox().hiVect();
        Real* dat     = mf[mfmfi].dataPtr(comp);
        Real* rad     = &radius[level][mfmfi.index()];

        FORT_RADMPY(dat,ARLIM(lo),ARLIM(hi),domlo,domhi,&ngrow,
                    rad,&nr,&bogus_value);
    }
}

//
// Divide by a radius for r-z coordinates.
//

void
Projection::radDiv (int       level,
                    MultiFab& mf,
                    int       comp)
{
    BL_ASSERT(comp >= 0 && comp < mf.nComp());
    BL_ASSERT(radius_grow >= mf.nGrow());

    int ngrow = mf.nGrow();
    int nr    = radius_grow;

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

    for (MFIter mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
    {
        BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

        const int* lo  = mfmfi.validbox().loVect();
        const int* hi  = mfmfi.validbox().hiVect();
        Real*      dat = mf[mfmfi].dataPtr(comp);
        Real*      rad = &radius[level][mfmfi.index()];

        FORT_RADDIV(dat,ARLIM(lo),ARLIM(hi),domlo,domhi,&ngrow,
                    rad,&nr,&bogus_value);
    }
}

//
// This projects the initial vorticity field (stored in pressure)
// to define an initial velocity field.
//

void
Projection::initialVorticityProject (int c_lev)
{
#if (BL_SPACEDIM == 2)
    int f_lev = finest_level;

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "initialVorticityProject: levels = "
                  << c_lev
                  << "  "
                  << f_lev << std::endl;
    }
    //
    // Set up projector bndry just for this projection.
    //
    const int*      lo_bc = phys_bc->lo();
    const int*      hi_bc = phys_bc->hi();
    const Geometry& geom  = parent->Geom(0);

    RegType proj_bc[BL_SPACEDIM][2];

    proj_bc[0][0] = outflow;
    proj_bc[0][1] = outflow;
    if (geom.isPeriodic(0))
    {
        proj_bc[0][0] = periodic;
        proj_bc[0][1] = periodic;
    }
    proj_bc[1][0] = outflow;
    proj_bc[1][1] = outflow;

    if (geom.isPeriodic(1))
    {
        proj_bc[1][0] = periodic;
        proj_bc[1][1] = periodic;
    }

    delete projector_bndry;

#ifdef BL_USE_HGPROJ_SERIAL
    projector_bndry = new inviscid_fluid_boundary_class(proj_bc);
#else
    projector_bndry = new inviscid_fluid_boundary(proj_bc);
#endif

    bldSyncProject();

    MultiFab* vel[MAX_LEV] = {0};

    PArray<MultiFab> p_real(f_lev+1,PArrayManage);
    PArray<MultiFab> s_real(f_lev+1,PArrayManage);

    for (int lev = c_lev; lev <= f_lev; lev++)
    {
        MultiFab& P_new  = LevelData[lev].get_new_data(Press_Type);
        const int nghost = 1;
        s_real.set(lev,new MultiFab(LevelData[lev].boxArray(),1,nghost));
        s_real[lev].setVal(1,nghost);
        p_real.set(lev,new MultiFab(P_new.boxArray(),1,nghost));
        p_real[lev].setVal(0,nghost);
    }
    //
    // Set up outflow bcs.
    //
    const Array<BoxArray>& full_mesh = sync_proj->mesh();

    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> rhs_real(f_lev+1,PArrayManage);

    for (int n = 0; n < BL_SPACEDIM; n++)
        u_real[n].resize(f_lev+1,PArrayManage);

    for (int lev = c_lev; lev <= f_lev; lev++)
    {
        for (int n = 0; n < BL_SPACEDIM; n++)
        {
            u_real[n].set(lev, new MultiFab(full_mesh[lev], 1, 1));
            u_real[n][lev].setVal(0);
        }
        //
        // The vorticity is stored in the new pressure variable for now.
        //
        MultiFab& P_new = LevelData[lev].get_new_data(Press_Type);

        rhs_real.set(lev, new MultiFab(P_new.boxArray(), 1, 1));

        for (MFIter mfi(rhs_real[lev]); mfi.isValid(); ++mfi)
        {
            rhs_real[lev][mfi].setVal(0);
            rhs_real[lev][mfi].copy(P_new[mfi], 0, 0);
        }
    }

    MultiFab* sync_resid_crse = 0;
    MultiFab* sync_resid_fine = 0;
    //
    // Project.
    //
    const Real* dx_lev = parent->Geom(f_lev).CellSize();
    const bool  use_u  = false;
    sync_proj->manual_project(u_real,p_real,rhs_real,null_amr_real,s_real,
                              sync_resid_crse, sync_resid_fine, 
                              parent->Geom(c_lev), 
                              use_u,(Real*)dx_lev,
                              proj_tol,c_lev,f_lev,proj_abs_tol);

    const int idx[2] = {1, 0};

    for (int lev = c_lev; lev <= f_lev; lev++)
    {
        vel[lev] = &LevelData[lev].get_new_data(State_Type);
        //
        // Note: Here u_real from projection is -grad(phi), but if
        //  phi is the stream function, u=dphi/dy, v=-dphi/dx
        //
        (u_real[Yvel][lev]).mult(-1,0,1);

        for (int n = 0; n < BL_SPACEDIM; n++)
        {
            for (MFIter mfi(*vel[lev]); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.validbox();
                if (add_vort_proj)
                {
                    (*vel[lev])[mfi].plus(u_real[n][lev][mfi],box,0,Xvel+idx[n], 1);
                }
                else
                {
                    (*vel[lev])[mfi].copy(u_real[n][lev][mfi],box,0,box,Xvel+idx[n], 1);
                }
            }
        }
    }
    //
    // Reset the boundary conditions for all the other projections.
    //
    setUpBcs();
    //
    // Remove the sync projector built with these bc's. 
    //
    delete sync_proj;

    sync_proj = 0;
#else
    BoxLib::Error("Projection::initialVorticityProject(): not implented yet for 3D");
#endif
}

void 
Projection::putDown (MultiFab**         phi,
                     MultiFab&          phi_fine_strip,
                     int                c_lev,
                     int                f_lev,
                     const Orientation* outFaces,
                     int                numOutFlowFaces,
                     int                ncStripWidth)
{
    //
    // Put down to coarser levels.
    //
    const int nCompPhi = phi_fine_strip.nComp();
    const int nGrow    = phi_fine_strip.nGrow();
    IntVect ratio      = IntVect::TheUnitVector();

    for (int lev = f_lev-1; lev >= c_lev; lev--) {

      ratio *= parent->refRatio(lev);
      const Box& domainC = parent->Geom(lev).Domain();
      BoxList phiC_strip_bl(IndexType::TheNodeType());
      
      for (int iface = 0; iface < numOutFlowFaces; iface++) 
      {
        Box phiC_strip = BoxLib::surroundingNodes(BoxLib::bdryNode(domainC,outFaces[iface],ncStripWidth));
        phiC_strip_bl.push_back(phiC_strip);
      }
      BoxArray phiC_strip_ba(phiC_strip_bl);
      MultiFab phi_crse_strip(phiC_strip_ba, nCompPhi, nGrow);
      phi_crse_strip.setVal(0);
      
        for (MFIter finemfi(phi_fine_strip); finemfi.isValid(); ++finemfi)
        {
            const int i = finemfi.index();
            Box ovlp = BoxLib::coarsen(finemfi.validbox(),ratio) & phi_crse_strip.box(i);
            FORT_PUTDOWN (phi_crse_strip[i].dataPtr(),
                          ARLIM(phi_crse_strip[i].loVect()),
                          ARLIM(phi_crse_strip[i].hiVect()),
                          phi_fine_strip[i].dataPtr(),
                          ARLIM(phi_fine_strip[i].loVect()),
                          ARLIM(phi_fine_strip[i].hiVect()),
                          ovlp.loVect(),ovlp.hiVect(),ratio.getVect());
        }
        phi[lev]->copy(phi_crse_strip);
    }
}

void
Projection::getStreamFunction (PArray<MultiFab>& phi)
{
#if (BL_SPACEDIM == 2)
    int c_lev = 0;
    int f_lev = finest_level;

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "getStreamFunction: levels = "
                  << c_lev
                  << "  "
                  << f_lev << std::endl;
    }
    //
    // Set up projector bndry just for this projection.
    //
    const int*      lo_bc = phys_bc->lo();
    const int*      hi_bc = phys_bc->hi();
    const Geometry& geom  = parent->Geom(0);

    RegType proj_bc[BL_SPACEDIM][2];

    proj_bc[0][0] = outflow;
    proj_bc[0][1] = outflow;
    if (geom.isPeriodic(0))
    {
        proj_bc[0][0] = periodic;
        proj_bc[0][1] = periodic;
    }
    proj_bc[1][0] = outflow;
    proj_bc[1][1] = outflow;

    if (geom.isPeriodic(1))
    {
        proj_bc[1][0] = periodic;
        proj_bc[1][1] = periodic;
    }

    delete projector_bndry;

#ifdef BL_USE_HGPROJ_SERIAL
    projector_bndry = new inviscid_fluid_boundary_class(proj_bc);
#else
    projector_bndry = new inviscid_fluid_boundary(proj_bc);
#endif

    bldSyncProject();

    MultiFab* vel[MAX_LEV] = {0};

    PArray<MultiFab> p_real(f_lev+1,PArrayManage);
    PArray<MultiFab> s_real(f_lev+1,PArrayManage);

    const int nghost = 1;
    for (int lev = c_lev; lev <= f_lev; lev++)
    {
        s_real.set(lev,new MultiFab(LevelData[lev].boxArray(),1,nghost));
        s_real[lev].setVal(1,nghost);
        p_real.set(lev,&phi[lev]);
    }
    //
    // Set up outflow bcs.
    //
    const Array<BoxArray>& full_mesh = sync_proj->mesh();

    PArray<MultiFab> u_real[BL_SPACEDIM];

    for (int n = 0; n < BL_SPACEDIM; n++)
        u_real[n].resize(f_lev+1,PArrayManage);
    //
    // Copy the velocity field into u_real.
    //
    for (int lev = c_lev; lev <= f_lev; lev++) 
    {
      vel[lev] = &LevelData[lev].get_new_data(State_Type);
      for (int n = 0; n < BL_SPACEDIM; n++) 
      {
        u_real[n].set(lev, new MultiFab(full_mesh[lev], 1, 1));
        for (int i = 0; i < full_mesh[lev].size(); i++) 
        {
          vel[lev] = &LevelData[lev].get_new_data(State_Type);
          u_real[n][lev][i].copy((*vel[lev])[i], n, 0);
        }
      }
    }
    //
    // Project.
    //
    const Real* dx_lev = parent->Geom(f_lev).CellSize();
    sync_proj->stream_func_project(u_real,p_real,s_real,
                                   (Real*)dx_lev, proj_tol,
                                   c_lev,f_lev,proj_abs_tol);

    //
    // Reset the boundary conditions for all the other projections.
    //
    setUpBcs();

    //
    // Remove the sync projector built with these bc's. 
    //
    delete sync_proj;
    sync_proj = 0;

#else
    BoxLib::Error("Projection::getStreamFunction(): not implented yet for 3D");
#endif
}

//
// Given a nodal pressure P compute the pressure gradient at the
// contained cell centers.

void
Projection::getGradP (FArrayBox& p_fab,
                      FArrayBox& gp,
                      const Box& gpbox_to_fill,
                      const Real* dx)
{
    //
    // Test to see if p_fab contains gpbox_to_fill
    //
    BL_ASSERT(BoxLib::enclosedCells(p_fab.box()).contains(gpbox_to_fill));

    const int*  plo    = p_fab.loVect();
    const int*  phi    = p_fab.hiVect();
    const int*  glo    = gp.box().loVect();
    const int*  ghi    = gp.box().hiVect();
    const int*   lo    = gpbox_to_fill.loVect();
    const int*   hi    = gpbox_to_fill.hiVect();
    const Real* p_dat  = p_fab.dataPtr();
    const Real* gp_dat = gp.dataPtr();

#if (BL_SPACEDIM == 2)
    int is_full = 0;
    if (hg_stencil == holy_grail_amr_multigrid::full)  is_full = 1;
    FORT_GRADP(p_dat,ARLIM(plo),ARLIM(phi),gp_dat,ARLIM(glo),ARLIM(ghi),lo,hi,dx,
               &is_full);
#elif (BL_SPACEDIM == 3)
    FORT_GRADP(p_dat,ARLIM(plo),ARLIM(phi),gp_dat,ARLIM(glo),ARLIM(ghi),lo,hi,dx);
#endif
}

void
Projection::set_outflow_bcs (int        which_call,
                             MultiFab** phi, 
                             MultiFab** Vel_in,
                             MultiFab** Divu_in,
                             MultiFab** Sig_in,
                             int        c_lev,
                             int        f_lev,
                             int        have_divu)
{
    BL_ASSERT((which_call == INITIAL_VEL  ) || 
              (which_call == INITIAL_PRESS) || 
              (which_call == INITIAL_SYNC ) ||
              (which_call == LEVEL_PROJ   ) );

    if (which_call != LEVEL_PROJ)
      BL_ASSERT(c_lev == 0);

    std::cout << "...setting outflow bcs for the nodal projection ... " << std::endl;

    bool        hasOutFlow;
    Orientation outFaces[2*BL_SPACEDIM];
    Orientation outFacesAtThisLevel[MAXLEV][2*BL_SPACEDIM];

    int fine_level[2*BL_SPACEDIM];

    int numOutFlowFacesAtAllLevels;
    int numOutFlowFaces[MAXLEV];
    OutFlowBC::GetOutFlowFaces(hasOutFlow,outFaces,phys_bc,numOutFlowFacesAtAllLevels);

    //
    // Get 2-wide cc box, state_strip, along interior of top. 
    // Get 1-wide nc box, phi_strip  , along top.
    //
    const int ccStripWidth = 2;

    const int nCompPhi    = 1;
    const int srcCompVel  = Xvel;
    const int srcCompDivu = 0;
    const int   nCompVel  = BL_SPACEDIM;
    const int   nCompDivu = 1;

    //
    // Determine the finest level such that the entire outflow face is covered
    // by boxes at this level (skip if doesnt touch, and bomb if only partially
    // covered).
    //
    Box state_strip[MAXLEV][2*BL_SPACEDIM];

    int icount[MAXLEV];
    for (int i=0; i < MAXLEV; i++) icount[i] = 0;

    //
    // This loop is only to define the number of outflow faces at each level.
    //
    Box temp_state_strip;
    for (int iface = 0; iface < numOutFlowFacesAtAllLevels; iface++) 
    {
      const int outDir    = outFaces[iface].coordDir();

      fine_level[iface] = -1;
      for (int lev = f_lev; lev >= c_lev; lev--)
      {
        Box domain = parent->Geom(lev).Domain();

        if (outFaces[iface].faceDir() == Orientation::high)
        {
            temp_state_strip
                = Box(BoxLib::adjCellHi(domain,outDir,ccStripWidth)).shift(outDir,-ccStripWidth);
        }
        else
        {
            temp_state_strip
                = Box(BoxLib::adjCellLo(domain,outDir,ccStripWidth)).shift(outDir,ccStripWidth);
        }
        // Grow the box by one tangentially in order to get velocity bc's.
        for (int dir = 0; dir < BL_SPACEDIM; dir++) 
          if (dir != outDir) temp_state_strip.grow(dir,1);

        const BoxArray& Lgrids               = parent->getLevel(lev).boxArray();
        const Box       valid_state_strip    = temp_state_strip & domain;
        const BoxArray  uncovered_outflow_ba = BoxLib::complementIn(valid_state_strip,Lgrids);

        BL_ASSERT( !(uncovered_outflow_ba.size() &&
                     BoxLib::intersect(Lgrids,valid_state_strip).size()) );

        if ( !(uncovered_outflow_ba.size()) && fine_level[iface] == -1) {
            int ii = icount[lev];
            outFacesAtThisLevel[lev][ii] = outFaces[iface];
            state_strip[lev][ii] = temp_state_strip;
            fine_level[iface] = lev;
            icount[lev]++;
        }

      // end level loop
      }
    // end iface loop
    }

    for (int lev = f_lev; lev >= c_lev; lev--) {
      numOutFlowFaces[lev] = icount[lev];
    }

    NavierStokes* ns0 = dynamic_cast<NavierStokes*>(&LevelData[c_lev]);
    BL_ASSERT(!(ns0 == 0));
   
    int Divu_Type, Divu;
    Real gravity;

    if (which_call == INITIAL_SYNC || which_call == INITIAL_VEL)
    {
      gravity = 0.;
      if (!LevelData[c_lev].isStateVariable("divu", Divu_Type, Divu))
        BoxLib::Error("Projection::set_outflow_bcs: No divu.");
    }

    if (which_call == INITIAL_PRESS || which_call == LEVEL_PROJ)
    {
      gravity = ns0->getGravity();
      if (!LevelData[c_lev].isStateVariable("divu", Divu_Type, Divu) &&
          (gravity == 0.) )
        BoxLib::Error("Projection::set_outflow_bcs: No divu or gravity.");
    }

    for (int lev = c_lev; lev <= f_lev; lev++) 
    {
      if (numOutFlowFaces[lev] > 0) 
        set_outflow_bcs_at_level (which_call,lev,c_lev,
                                  state_strip[lev],
                                  outFacesAtThisLevel[lev],
                                  numOutFlowFaces[lev],
                                  phi,
                                  Vel_in[lev],
                                  Divu_in[lev],
                                  Sig_in[lev],
                                  have_divu,
                                  gravity);
                                  
    }

}

void
Projection::set_outflow_bcs_at_level (int          which_call,
                                      int          lev,
                                      int          c_lev,
                                      Box*         state_strip,
                                      Orientation* outFacesAtThisLevel,
                                      int          numOutFlowFaces,
                                      MultiFab**   phi, 
                                      MultiFab*    Vel_in,
                                      MultiFab*    Divu_in,
                                      MultiFab*    Sig_in,
                                      int          have_divu,
                                      Real         gravity)
{
    NavierStokes* ns = dynamic_cast<NavierStokes*>(&LevelData[lev]);
    BL_ASSERT(!(ns == 0));

    Box domain = parent->Geom(lev).Domain();

    BoxList   phi_strip_bl(IndexType::TheNodeType());
    BoxList state_strip_bl;

    const int nGrow        = 0;
    const int nCompPhi     = 1;
    const int ncStripWidth = 1;

    std::cout << "NUM FACES " << numOutFlowFaces << std::endl;

    for (int iface = 0; iface < numOutFlowFaces; iface++) {
        Box phi_strip = 
         BoxLib::surroundingNodes(BoxLib::bdryNode(domain,
                   outFacesAtThisLevel[iface],ncStripWidth));
        phi_strip_bl.push_back(phi_strip);
        state_strip_bl.push_back(state_strip[iface]);
    }

    BoxArray   phi_strip_ba(  phi_strip_bl);
    BoxArray state_strip_ba(state_strip_bl);

    MultiFab phi_fine_strip(phi_strip_ba,nCompPhi,nGrow);
    phi_fine_strip.setVal(0);

    MultiFab  rho(state_strip_ba,          1,nGrow);
    MultiFab dsdt(state_strip_ba,          1,nGrow);
    MultiFab dudt(state_strip_ba,BL_SPACEDIM,nGrow);

     rho.setVal(1.e200);
    dudt.setVal(1.e200);
    dsdt.setVal(1.e200);

    for (MFIter mfi(rho); mfi.isValid(); ++mfi)
       Sig_in->copy(rho[mfi.index()]);

    ProjOutFlowBC projBC;
    if (which_call == INITIAL_PRESS) {

      projBC.computeRhoG(rho,phi_fine_strip,
                         parent->Geom(lev),
                         outFacesAtThisLevel,numOutFlowFaces,gravity);

    } else {

      if (have_divu) {
        Vel_in->FillBoundary();

//      Build a new MultiFab for which the cells outside the domain
//        are in the valid region instead of being ghost cells, so that
//        we can copy these values into the dudt array.
        BoxList grown_vel_bl;
        for (int i = 0; i < Vel_in->size(); i++)
           grown_vel_bl.push_back(BoxLib::grow(Vel_in->boxArray()[i],1));
        BoxArray grown_vel_ba(grown_vel_bl);
        MultiFab grown_vel(grown_vel_ba,BL_SPACEDIM,0);
        for (MFIter vmfi(*Vel_in); vmfi.isValid(); ++vmfi)
           grown_vel[vmfi.index()].copy((*Vel_in)[vmfi.index()]);
        
        for (MFIter mfi(dudt); mfi.isValid(); ++mfi)
        {
          grown_vel.copy(dudt[mfi.index()]);
           Divu_in->copy(dsdt[mfi.index()]);
        }

      } else {
        for (MFIter mfi(dudt); mfi.isValid(); ++mfi)
        {
          Vel_in->copy(dudt[mfi.index()]);
          dsdt[mfi.index()].setVal(0.);
        }
      }


      projBC.computeBC(&dudt,dsdt,rho,phi_fine_strip,
                     parent->Geom(lev),
                     outFacesAtThisLevel,
                     numOutFlowFaces,gravity);

    }

    phi[lev]->copy(phi_fine_strip);

    if (lev > c_lev) 
      putDown(phi,phi_fine_strip,c_lev,lev,outFacesAtThisLevel,
              numOutFlowFaces,ncStripWidth);
}
