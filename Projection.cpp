//BL_COPYRIGHT_NOTICE

//
// $Id: Projection.cpp,v 1.88 1999-06-03 18:17:12 almgren Exp $
//

#ifdef BL_T3E
#include <List.H>
#endif

#include <Misc.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <MacOpProjDrivers.H>
#include <NavierStokes.H>
#include <Projection.H>
#include <RunStats.H>
#include <SYNCREG_F.H>
#include <PROJECTION_F.H>
#include <NAVIERSTOKES_F.H>
#include <hg_projector.H>

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

//
// This is a temporary function.
// Eventually this will be moved to a boundary condition class.
//

static
void
getOutFlowFace (bool&        haveOutFlow,
                Orientation& outFace,
                BCRec*       _phys_bc)
{
    haveOutFlow = false;

    int numOutFlowBC = 0;

    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
        if (_phys_bc->lo(idir) == Outflow)
        {
            haveOutFlow = true;
            outFace = Orientation(idir,Orientation::low);
            numOutFlowBC++;
        }
        if (_phys_bc->hi(idir) == Outflow)
        {
            haveOutFlow = true;
            outFace = Orientation(idir,Orientation::high);
            numOutFlowBC++;
        }
    }

    if (numOutFlowBC > 1)
        BoxLib::Error("currently only allowed one outflow bc");
}

//
// This is a temporary function.
// Eventually this will be moved to a boundary condition class.
//

static
bool
hasOutFlowBC (BCRec* _phys_bc)
{
  bool has_out_flow = false;
  
  int numOutFlowBC = 0;

  for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
      if (_phys_bc->lo(idir) == Outflow)
        {
	  has_out_flow = true;
	  numOutFlowBC++;
        }
      if (_phys_bc->hi(idir) == Outflow)
        {
	  has_out_flow = true;
	  numOutFlowBC++;
        }
    }
  
  if (numOutFlowBC > 1) 
    BoxLib::Error("currently only allowed one outflow bc");
  
  return has_out_flow;
}

// this is a temporary function. it will go away when the bc are generalized.
#if (BL_SPACEDIM != 2)
static 
Real
checkDivU(Amr* parent,
	  PArray<AmrLevel>& LevelData,
	  BCRec* phys_bc,
	  int f_lev,
	  int c_lev,
	  Real cur_divu_time)
{
    // Get 3-wide cc box, state_strip, along top,
    bool hasOutFlow;
    Orientation _outFace;
    getOutFlowFace(hasOutFlow,_outFace,phys_bc);

    const Box& domain      = parent->Geom(f_lev).Domain();
    const int outDir       = _outFace.faceDir();
    const int ccStripWidth = 3;
    Box state_strip;
    if (_outFace.isLow()) {
      state_strip = Box(::adjCellLo(domain,outDir,ccStripWidth)).shift(outDir,ccStripWidth);
    } else {
      state_strip = Box(::adjCellHi(domain,outDir,ccStripWidth)).shift(outDir,-ccStripWidth);
    }
    const int nGrow       = 0;
    const int srcCompDivu = 0,      nCompDivu = 1;

    BoxArray state_strip_ba(&state_strip,1);
    MultiFab cc_MultiFab(state_strip_ba, 1, nGrow, Fab_noallocate);

    int Divu_Type, Divu;
    if (!LevelData[c_lev].isStateVariable("divu", Divu_Type, Divu)) 
        BoxLib::Error("Projection::set_initial_projection_outflow_bcs(): Divu not found");
    
    FillPatchIterator divuFpi(LevelData[f_lev],cc_MultiFab,nGrow,
                          cur_divu_time,Divu_Type,srcCompDivu,nCompDivu);
    FArrayBox divu_fab(state_strip,1);
    for ( ; divuFpi.isValid(); ++divuFpi)
    {
      divu_fab.copy(divuFpi());
    }
    REAL norm_divu = divu_fab.norm(1,srcCompDivu,nCompDivu);
    return norm_divu;
}
#endif

#define BogusValue 1.e20
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

static RegType project_bc [] =
{
    interior, inflow, outflow, refWall, refWall, refWall
};

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
#ifdef BL_T3E
    //
    // Force instantiation of List<int>::clear() for the T3E.
    //
    List<int> tempIntList;
    tempIntList.clear();
#endif

    read_params();

    if (verbose && ParallelDescriptor::IOProcessor()) 
        cout << "Creating projector\n";

#if BL_SPACEDIM == 2
    if (CoordSys::IsRZ())
        amr_multigrid::SetRZ();
#endif
    setUpBcs();
    sync_proj = 0;
}

Projection::~Projection ()
{
    if (verbose && ParallelDescriptor::IOProcessor()) 
    {
        cout << "Deleting projector\n";
    }
    delete sync_proj;
    delete projector_bndry;
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

    pp.query("divu_minus_s_factor",divu_minus_s_factor);

    pp.query("rho_wgt_vel_proj",rho_wgt_vel_proj);

    pp.query("do_outflow_bcs",do_outflow_bcs);
}

void 
Projection::setUpBcs ()
{
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

    projector_bndry = new inviscid_fluid_boundary_class(proj_bc);
}

//
// Install a level of the projection.
//

void
Projection::install_level (int           level,
                           AmrLevel*     level_data,
                           PArray<Real>* _radius )
{
    if (verbose && ParallelDescriptor::IOProcessor()) 
        cout << "Installing projector level " << level << '\n';

    finest_level = parent->finestLevel();

    if (level > LevelData.length() - 1) 
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
        cout << "bldSyncProject:: amr_mesh = \n";
        amr_multigrid::mesh_write(amesh, gen_ratio, fdomain, cout);
    }

#ifdef BL_USE_HGPROJ_SERIAL
    sync_proj = new holy_grail_amr_projector(amesh, gen_ratio, fdomain,
                                             0, finest_level, finest_level,
                                             *projector_bndry, P_code);
#else
    sync_proj = new holy_grail_amr_projector(amesh, gen_ratio, fdomain,
                                             0, finest_level, finest_level,
                                             *projector_bndry,
					     holy_grail_amr_multigrid::cross,
                                             P_code);
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
    if (ParallelDescriptor::IOProcessor() && verbose)
	cout << "... level projector at level " << level << '\n';
    
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
    // Compute Ustar/dt as input to projector for proj_0,
    //         Ustar/dt + Gp                  for proj_2,
    //         (Ustar-Un)/dt for not proj_0 or proj_2 (ie the original).
    //
    // Compute DU/dt for proj_0,
    //         DU/dt for proj_2,
    //         (DU-DU_old)/dt for not proj_0 or proj_2 (ie the original).
    //

    MultiFab *divusource, *divuold;
    divusource = divuold = 0;
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
        for (MultiFabIterator U_newmfi(U_new); U_newmfi.isValid(); ++U_newmfi) 
        {
            DependentMultiFabIterator U_oldmfi(U_newmfi, U_old);
            ConvertUnew(U_newmfi(),U_oldmfi(),dt,U_newmfi.validbox());
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

    REAL prev_pres_time = cur_pres_time - dt;
    FArrayBox Gp;

    if (proj_2)
    {
      FillPatchIterator P_oldfpi(LevelData[level],P_old,1,
                                 prev_pres_time,Press_Type,0,1);
      int n_ghost = 1;
      rho_half->invert(1.0,n_ghost);

      for ( ; P_oldfpi.isValid(); ++P_oldfpi)
      {
        ns->getGradP(P_oldfpi(),Gp,grids[P_oldfpi.index()],n_ghost);

        Gp.mult((*rho_half)[P_oldfpi.index()],0,0,1);
        Gp.mult((*rho_half)[P_oldfpi.index()],0,1,1);
#if (BL_SPACEDIM==3)
        Gp.mult((*rho_half)[P_oldfpi.index()],0,2,1);
#endif

        DependentMultiFabIterator U_newmfi(P_oldfpi, U_new);
        U_newmfi().plus(Gp,0,0,BL_SPACEDIM);
      }

      rho_half->invert(1.0,n_ghost);
    }

    //
    // Set boundary values for P_new, to increment, if applicable
    //
    if (level != 0)
    {
	LevelData[level].FillCoarsePatch(P_new,0,cur_pres_time,Press_Type,0,1);
        if (!proj_2) 
          P_new.minus(P_old,0,1,0); // Care about nodes on box boundary
    }
    const int nGrow = (level == 0  ?  0  :  -1);
    for (MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi)
    {
        P_newmfi().setVal(0.0,::grow(P_newmfi.validbox(),nGrow),0,1);
        // Also, zero fine-fine nodes?
    }
    //
    // Overwrite IC with outflow Dirichlet, if applicable
    //
    if(hasOutFlowBC(phys_bc) && have_divu && do_outflow_bcs) 
    {
	set_level_projector_outflow_bcs(level,time+0.5*dt,P_new,U_new,*divusource);
    }
    //
    // Scale the projection variables.
    //
    rho_half->setBndry(BogusValue);
    scaleVar(rho_half, 1, &U_new, grids, level);
    if (have_divu)
        radMult(level,*divusource,0);
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
    int rz_flag = (CoordSys::IsRZ() ? 1 : 0);

    if (do_sync_proj) 
    {
        int isrz   = CoordSys::IsRZ();
        int bcxlo  = phys_bc->lo(0);
        int bcxhi  = phys_bc->hi(0);
        int lowfix = (isrz==1 && bcxlo!=FOEXTRAP && bcxlo!=HOEXTRAP);
        int hifix  = (isrz==1 && bcxhi!=FOEXTRAP && bcxhi!=HOEXTRAP);

        if (level < finest_level) 
        {
            //
            // Init sync registers between level and level+1.
            //
            crse_sync_reg->CrseDVInit(U_new,geom,rz_flag,bc);
            if (have_divu) 
                crse_sync_reg->CrseDsdtAdd(*divusource,geom,rz_flag,bc,lowfix,hifix);
        } 
        if (level > 0 && (((proj_0 || proj_2) &&
                           iteration == crse_dt_ratio) ||
                          (!proj_0 && !proj_2)))
        {
            //
            // Increment sync registers between level and level-1.
            //
            const Real invrat         = 1.0/(double)crse_dt_ratio;
            const Geometry& crse_geom = parent->Geom(level-1);
            fine_sync_reg->FineDVAdd(U_new,dx,crse_geom,rz_flag,bc,invrat);
            if (have_divu) 
                fine_sync_reg->FineDsdtAdd(*divusource,geom,crse_geom,rz_flag,
                                           bc,lowfix,hifix,invrat);
        }
    }
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
        u_real[n].resize(level+1);
        u_real[n].set(level, new MultiFab(grids, 1, 1));
        for (MultiFabIterator u_realmfi(u_real[n][level]); u_realmfi.isValid();
            ++u_realmfi)
        {
            DependentMultiFabIterator U_newmfi(u_realmfi, U_new);
            u_realmfi().copy(U_newmfi(), n, 0);
        }
    }
    p_real.set(level, &P_new);
    s_real.set(level, rho_half);
    //
    // Project
    //
    if (!have_divu) 
    {
        sync_proj->project(u_real, p_real, null_amr_real, s_real, (Real*)dx,
                           proj_tol, level, level, proj_abs_tol);
    } 
    else 
    {
        bool      use_u  = true;
        const int nghost = 1;
        divusource->mult(-1.0,0,1,nghost);

        PArray<MultiFab> rhs_real(level+1);
        rhs_real.set(level, divusource);
        sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real,
                                  s_real, use_u, (Real*)dx,
                                  proj_tol, level, level, proj_abs_tol);
    }
    //
    // Copy and delete u_real.
    //
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        for (MultiFabIterator u_realmfi(u_real[n][level]); u_realmfi.isValid();
            ++u_realmfi)
        {
            DependentMultiFabIterator U_newmfi(u_realmfi, U_new);
            U_newmfi().copy(u_realmfi(), 0, n);
        }
        delete u_real[n].remove(level);
    }
    //
    // Increment sync registers by adding contribution from the
    // full D(sig Gphi)
    //
    if (do_sync_proj) 
    {
        if (level < finest_level) 
        {
            //
            // Init sync registers between level and level+1.
            //
            crse_sync_reg->CrseLPhiAdd(P_new,*rho_half,geom,rz_flag);
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
            fine_sync_reg->FineLPhiAdd(P_new,*rho_half,dx,crse_geom,rz_flag,invrat);
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

    delete divusource;
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
        cout << "... filterP at level " << level << endl;

    int             rzflag   = CoordSys::IsRZ();
    const Real*     dx       = geom.CellSize();
    Box             ndomain  = ::surroundingNodes(geom.Domain());
    const int*      ndlo     = ndomain.loVect();
    const int*      ndhi     = ndomain.hiVect();
    const BoxArray& grids    = LevelData[level].boxArray();
    const BoxArray& P_grids  = P_old.boxArray();
    const int*      phys_lo  = phys_bc->lo();
    const int*      phys_hi  = phys_bc->hi();
    MultiFab*       rhs_cc   = new MultiFab(P_grids,1,1,Fab_allocate);
    MultiFab*       temp_phi = new MultiFab(P_grids,1,1,Fab_allocate);
    MultiFab*       temp_rho = new MultiFab(grids,1,1,Fab_allocate);
    MultiFab*       temp_vel = new MultiFab(grids,BL_SPACEDIM,1,Fab_allocate);

    MultiFab* divuold = 0;

    BL_ASSERT(grids.length() == P_grids.length());
    
    temp_phi->setVal(0);
    temp_rho->setVal(0);
    rhs_cc->setVal(0);
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

    for (MultiFabIterator mfi(*temp_rho); mfi.isValid(); ++mfi)
    {
        const int  k    = mfi.index();
        FArrayBox& sfab = (*temp_rho)[k];
        const int* s_lo = sfab.loVect();
        const int* s_hi = sfab.hiVect();
        FArrayBox& pfab = (*temp_phi)[k];
        const int* p_lo = pfab.loVect();
        const int* p_hi = pfab.hiVect();
        const int* r_lo = (*rhs_cc)[k].loVect();
        const int* r_hi = (*rhs_cc)[k].hiVect();
        const int* n_lo = P_grids[k].loVect();
        const int* n_hi = P_grids[k].hiVect();

        FORT_FILTRHS(pfab.dataPtr(),ARLIM(p_lo),ARLIM(p_hi),
                     sfab.dataPtr(),ARLIM(s_lo),ARLIM(s_hi),
                     (*rhs_cc)[k].dataPtr(),ARLIM(r_lo),ARLIM(r_hi),
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
                // two (only for ref-wall and inflow) and set to zero at
                // outflow.
                //
                Box blo = P_grids[k] & ndomlo;
                Box bhi = P_grids[k] & ndomhi;

                if (blo.ok())
                {
                    if (phys_lo[dir] == Outflow) 
                        (*rhs_cc)[k].setVal(0,blo,0,1);
                    else
                        (*rhs_cc)[k].mult(2,blo,0,1);
                }
                if (bhi.ok())
                {
                    if (phys_hi[dir] == Outflow) 
                        (*rhs_cc)[k].setVal(0,bhi,0,1);
                    else
                        (*rhs_cc)[k].mult(2,bhi,0,1);
                }
            } 
        }
    }

    if (have_divu)
    {
        NavierStokes* ns = dynamic_cast<NavierStokes*>(&parent->getLevel(level));

        BL_ASSERT(!(ns == 0));

        MultiFab* divuold = ns->getDivCond(1,time);

        for (MultiFabIterator mfi(*divuold); mfi.isValid(); ++mfi)
        {
            mfi().mult(1.0/dt,0,1);
        }

        const int nghost = 1;
        radMult(level,*divuold,0);
        rhs_cc->plus(*divuold,0,1,nghost);
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
        u_real[n].resize(level+1);
        u_real[n].set(level, new MultiFab(grids, 1, 1));
        for (MultiFabIterator u_realmfi(u_real[n][level]); u_realmfi.isValid();
            ++u_realmfi)
        {
            DependentMultiFabIterator velmfi(u_realmfi, *temp_vel);
            u_realmfi().copy(velmfi(), n, 0);
        }
    }
    p_real.set(level, temp_phi);
    s_real.set(level, temp_rho);
    rhs_real.set(level, rhs_cc);
    //
    // Project ...
    //
    const bool use_u = 1;
    sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                              use_u, (Real*)dx,
                              filter_factor, level, level, proj_abs_tol);

    for (int n = 0; n < BL_SPACEDIM; n++)
        delete u_real[n].remove(level);
    //
    // Reset state + pressure data
    //
    AddPhi(P_new, *temp_phi, grids);
    //
    // Unscale the projection variables.
    //
    rescaleVar(rho_half, 0, &U_old, grids, level);

    delete temp_phi;
    delete temp_rho;
    delete temp_vel;
    delete divuold;
    delete rhs_cc;
}

// 
//  This function is an attempt at improving the creation of fine
//  pressure data after regridding
//

void
Projection::harmonic_project (int             level,
                              Real            dt,
                              Real            cur_pres_time,
                              const Geometry& geom,
                              MultiFab&       P_old)
{
    BL_ASSERT(level != 0);

    if (verbose && ParallelDescriptor::IOProcessor()) 
        cout << "... harmonic projector\n";

    if (sync_proj == 0)
        bldSyncProject();

    const BoxArray& grids   = LevelData[level].boxArray();
    const BoxArray& P_grids = P_old.boxArray();

    MultiFab* rhs      = new MultiFab(P_grids,1,1,Fab_allocate);
    MultiFab* harm_phi = new MultiFab(P_grids,1,1,Fab_allocate);
    MultiFab* temp_phi = new MultiFab(P_grids,1,1,Fab_allocate);
    MultiFab* rho      = new MultiFab(grids,1,1,Fab_allocate);
    MultiFab* harm_vel = new MultiFab(grids,BL_SPACEDIM,1,Fab_allocate);

    rhs->setVal(0);
    harm_phi->setVal(0);
    harm_vel->setVal(0);
    rho->setVal(1);

    harm_phi->setBndry(BogusValue);

    const Real prev_pres_time = cur_pres_time - dt;

    LevelData[level].FillCoarsePatch(*temp_phi,0,prev_pres_time,Press_Type,0,1);
    LevelData[level].FillCoarsePatch(*harm_phi,0,cur_pres_time,Press_Type,0,1);

    for (MultiFabIterator phimfi(*temp_phi); phimfi.isValid(); ++phimfi)
    {
        DependentMultiFabIterator harm_phimfi(phimfi, *harm_phi);
        harm_phimfi().minus(phimfi());
        Box tempbox(harm_phimfi().box());
        tempbox.grow(-2);
        harm_phimfi().setVal(0,tempbox,0,1);
    }

    delete temp_phi;

    rho->setBndry(BogusValue);
    scaleVar(rho,1,0,grids,level);

    const Real* dx = geom.CellSize();
    //
    // Build alias lib structures.
    //
    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> p_real(level+1), s_real(level+1), rhs_real(level+1);
    int n;
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        u_real[n].resize(level+1);
        u_real[n].set(level, new MultiFab(grids, 1, 1));
        for (MultiFabIterator u_realmfi(u_real[n][level]); u_realmfi.isValid();
            ++u_realmfi)
        {
            DependentMultiFabIterator harm_velmfi(u_realmfi, *harm_vel);
            u_realmfi().copy(harm_velmfi(), n, 0);
        }
    }
    p_real.set(level, harm_phi);
    s_real.set(level, rho);
    rhs_real.set(level, rhs);
    //
    // Project
    //
    bool use_u = false;
    sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                              use_u, (Real*)dx,
                              proj_tol, level, level, proj_abs_tol);
    //
    // Copy and delete u_real.
    //
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        //
        // Copy unnecessary since harm_vel will be discarded anyway.
        //
        delete u_real[n].remove(level);
    }
    //
    // Unscale variables for harmonic projection.
    //
    rescaleVar(rho,1,0,grids,level);
    //
    // Update pressure.
    //
    AddPhi(P_old, *harm_phi, grids);

    delete rhs;
    delete harm_phi;
    delete harm_vel;
    delete rho;
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
                         int**           sync_bc,
                         const Geometry& geom,
                         const Real*     dx,
                         Real            dt_crse,
                         int             crse_dt_ratio)
{
    static RunStats stats("sync_project");

    stats.start();

    int rz_flag = (CoordSys::IsRZ() ? 1 : 0);

    if (verbose && ParallelDescriptor::IOProcessor()) 
    {
        cout << "SyncProject: level = "
             << c_lev
             << " correction to level "
             << finest_level << endl;
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
    // If this sync project is not at level 0 then we need to account for
    // the changes made here in the level c_lev velocity in the sync registers
    // going into the level (c_lev-1) sync project.
    //
    if (c_lev > 0) 
    {
        const Real      invrat    = 1.0/crse_dt_ratio;
        const Geometry& crsr_geom = parent->Geom(c_lev-1);
        crsr_sync_reg->CompDVAdd(*Vsync,sync_boxes,dx,geom,crsr_geom,
                                 rz_flag,sync_bc,invrat);
    }
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
        u_real[n].resize(c_lev+1);
        u_real[n].set(c_lev, new MultiFab(grids, 1, 1));
        for (MultiFabIterator u_realmfi(u_real[n][c_lev]); u_realmfi.isValid();
            ++u_realmfi)
        {
            DependentMultiFabIterator Vsyncmfi(u_realmfi, *Vsync);
            u_realmfi().copy(Vsyncmfi(), n, 0);
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
    bool use_u = true;
    sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                              use_u, (Real*)dx,
                              sync_tol, c_lev, c_lev, proj_abs_tol);
    //
    // Copy and delete u_real.
    //
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        for (MultiFabIterator u_realmfi(u_real[n][c_lev]); u_realmfi.isValid();
            ++u_realmfi)
        {
            DependentMultiFabIterator Vsyncmfi(u_realmfi, *Vsync);
            Vsyncmfi().copy(u_realmfi(), 0, n);
        }
        delete u_real[n].remove(c_lev);
    }
    //
    // If this sync project is not at level 0 then we need to account for
    // the changes made here in the level c_lev velocity in the sync registers
    // going into the level (c_lev-1) sync project.  Note that this must be
    // done before rho_half is scaled back.
    //
    if (c_lev > 0) 
    {
        const Real invrat = 1.0/(double)crse_dt_ratio;
        const Geometry& crsr_geom = parent->Geom(c_lev-1);
        crsr_sync_reg->CompLPhiAdd(phi,sig,sync_boxes,
                                   dx,geom,crsr_geom,rz_flag,invrat);
    }
    //
    //----------------- reset state + pressure data ---------------------
    //
    // Unscale the sync projection variables for rz.
    //
    rescaleVar(&sig,1,Vsync,grids,c_lev);
    //
    // Add projected Vsync to new velocity at this level & add phi to pressure.
    //
    AddPhi(pres, phi, grids);
    UpdateArg1(vel, dt_crse, *Vsync, BL_SPACEDIM, grids, 1);

    stats.end();
}

//
//  MULTI-LEVEL SYNC_PROJECT
//

void
Projection::MLsyncProject (int             c_lev,
                           MultiFab&       pres_crse,
                           MultiFab&       vel_crse,
                           MultiFab&       pres_fine,
                           MultiFab&       vel_fine,
                           MultiFab&       rho_crse,
                           MultiFab&       rho_fine,
                           MultiFab*       Vsync,
                           MultiFab&       V_corr,
                           MultiFab&       phi_fine,
                           SyncRegister*   rhs_sync_reg,
                           SyncRegister*   crsr_sync_reg,
                           int**           sync_bc,
                           const Real*     dx,
                           Real            dt_crse, 
                           IntVect&        ratio,
                           int             crse_dt_ratio,
                           const Geometry& fine_geom,
                           const Geometry& crse_geom,
                           int             have_divu)
{
    static RunStats stats("sync_project");

    stats.start();
    
    int lev;
    if (verbose && ParallelDescriptor::IOProcessor()) 
        cout << "SyncProject: levels = " << c_lev << ", " << c_lev+1 << '\n';
    
    int rz_flag = (CoordSys::IsRZ() ? 1 : 0);
    if (sync_proj == 0) bldSyncProject();
    //
    // Set up memory.
    //
    MultiFab *phi[MAXLEV];
    
    const BoxArray& grids      = LevelData[c_lev].boxArray();
    const BoxArray& fine_grids = LevelData[c_lev+1].boxArray();
    const BoxArray& Pgrids_crse = pres_crse.boxArray();

    phi[c_lev] = new MultiFab(Pgrids_crse,1,1,Fab_allocate);
    phi[c_lev]->setVal(0);
    
    MultiFab* crse_rhs = new MultiFab(Pgrids_crse,1,1,Fab_allocate);
    
    const BoxArray& Pgrids_fine = pres_fine.boxArray();
    phi[c_lev+1] = new MultiFab(Pgrids_fine,1,1,Fab_allocate);
    phi[c_lev+1]->setVal(0);

    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> p_real(c_lev+2), s_real(c_lev+2), crse_rhs_real(c_lev+1);

    BoxArray sync_boxes = pres_fine.boxArray();
    sync_boxes.coarsen(ratio);
    //
    // Set up RHS
    //
    rhs_sync_reg->InitRHS(*crse_rhs,crse_geom,phys_bc);
    
    Box P_finedomain(surroundingNodes(crse_geom.Domain()));
    P_finedomain.refine(ratio);
    if (Pgrids_fine[0] == P_finedomain) crse_rhs->setVal(0);
    
    //
    // Do necessary scaling
    //
    scaleVar(&rho_crse, 0, Vsync,   grids,      c_lev  );
    scaleVar(&rho_fine, 0, &V_corr, fine_grids, c_lev+1);
    //
    // If this sync project is not at level 0 then we need to account for
    // the changes made here in the level c_lev velocity in the sync registers
    // going into the level (c_lev-1) sync project.
    //
    if (c_lev > 0) 
    {
        const Real invrat = 1.0/(double)crse_dt_ratio;
        const Geometry& crsr_geom = parent->Geom(c_lev-1);
        crsr_sync_reg->CompDVAdd(*Vsync,sync_boxes,dx,
                                 crse_geom,crsr_geom,
                                 rz_flag,sync_bc,invrat);
    }
    //
    // Set up alias lib.
    //
    int n;
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        u_real[n].resize(c_lev+2);
        u_real[n].set(c_lev,   new MultiFab(grids,      1, 1));
        u_real[n].set(c_lev+1, new MultiFab(fine_grids, 1, 1));

        for (MultiFabIterator u_realmfi(u_real[n][c_lev]); u_realmfi.isValid();
            ++u_realmfi)
        {
            DependentMultiFabIterator Vsyncmfi(u_realmfi, *Vsync);
            u_realmfi().copy(Vsyncmfi(), n, 0);
        }
        for (MultiFabIterator u_realfinemfi(u_real[n][c_lev+1]); u_realfinemfi.isValid();
            ++u_realfinemfi)
        {
            DependentMultiFabIterator V_corrmfi(u_realfinemfi, V_corr);
            u_realfinemfi().copy(V_corrmfi(), n, 0);
        }

#ifdef BL_USE_HGPROJ_SERIAL
        restrict_level(u_real[n][c_lev], u_real[n][c_lev+1], ratio);
#else
        restrict_level(u_real[n][c_lev], u_real[n][c_lev+1], ratio, default_restrictor(), level_interface(), 0);
#endif
    }

    s_real.set(c_lev,   &rho_crse);
    s_real.set(c_lev+1, &rho_fine);

#ifdef BL_USE_HGPROJ_SERIAL
    restrict_level(s_real[c_lev], s_real[c_lev+1], ratio);
#else
    restrict_level(s_real[c_lev], s_real[c_lev+1], ratio, default_restrictor(), level_interface(), 0);
#endif

    p_real.set(c_lev,   phi[c_lev]);
    p_real.set(c_lev+1, phi[c_lev+1]);
    //
    // Note that crse_rhs_real is only built on the coarsest level.
    //
    crse_rhs_real.set(c_lev, crse_rhs);
    //
    // The Multilevel Projection
    // if use_u = 0, then solves DGphi = RHS
    // if use_u = 1, then solves DGphi = RHS + DV
    // both return phi and (V-Gphi) as V
    //
    const bool use_u     = true;
    const Real* dx_fine  = parent->Geom(c_lev+1).CellSize();

    sync_proj->manual_project(u_real, p_real, null_amr_real,
                              crse_rhs_real, s_real,
                              use_u, (Real*)dx_fine,
                              sync_tol, c_lev, c_lev+1, proj_abs_tol);
    //
    // Copy and delete u_real.
    //
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        for (MultiFabIterator u_realmfi(u_real[n][c_lev]); u_realmfi.isValid();
            ++u_realmfi)
        {
            DependentMultiFabIterator Vsyncmfi(u_realmfi, *Vsync);
            Vsyncmfi().copy(u_realmfi(), 0, n);
        }
        for (MultiFabIterator u_realfinemfi(u_real[n][c_lev+1]); u_realfinemfi.isValid();
            ++u_realfinemfi)
        {
            DependentMultiFabIterator V_corrmfi(u_realfinemfi, V_corr);
            V_corrmfi().copy(u_realfinemfi(), 0, n);
        }
        delete u_real[n].remove(c_lev);
        delete u_real[n].remove(c_lev+1);
    }
    //
    // If this sync project is not at levels 0-1 then we need to account for
    // the changes made here in the level c_lev velocity in the sync registers
    // going into the level (c_lev-1) sync project.  Note that this must be
    // done before rho_half is scaled back.
    //
    if (c_lev > 0) 
    {
        const Real invrat   = 1.0/(double)crse_dt_ratio;
        MultiFab*  phi_crse = phi[c_lev];
        const Geometry& crsr_geom = parent->Geom(c_lev-1);
        crsr_sync_reg->CompLPhiAdd(*phi_crse,rho_crse,sync_boxes,
                                   dx,crse_geom,crsr_geom,rz_flag,invrat);
    }
    //
    // Do necessary un-scaling.
    //
    rescaleVar(&rho_crse, 0, Vsync,   grids,      c_lev  );
    rescaleVar(&rho_fine, 0, &V_corr, fine_grids, c_lev+1);
    //
    // Add projected vel to new velocity and add phi to pressure .
    //
    AddPhi(pres_crse, *phi[c_lev],   grids     );
    AddPhi(pres_fine, *phi[c_lev+1], fine_grids);

    UpdateArg1(vel_crse, dt_crse, *Vsync, BL_SPACEDIM, grids,      1);
    UpdateArg1(vel_fine, dt_crse, V_corr, BL_SPACEDIM, fine_grids, 1);

    for (MultiFabIterator phimfi(*phi[c_lev+1]); phimfi.isValid(); ++phimfi) 
    {
        DependentMultiFabIterator phi_finemfi(phimfi, phi_fine);
        phi_finemfi().copy(phimfi(),0,0,1);
    }
    for (lev = c_lev; lev <= c_lev+1; lev++) 
        delete phi[lev];

    delete crse_rhs;

    stats.end();
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
    static RunStats stats("sync_project");

    stats.start();

    int lev;
    int f_lev = finest_level;
    if (verbose && ParallelDescriptor::IOProcessor()) 
    {
        cout << "initialVelocityProject: levels = " << c_lev
             << "  " << f_lev << '\n';
        if (rho_wgt_vel_proj) 
        {
            cout << "RHO WEIGHTED INITIAL VELOCITY PROJECTION\n";
        } 
        else 
        {
            cout << "CONSTANT DENSITY INITIAL VELOCITY PROJECTION\n";
        }
    }

    if (sync_proj == 0)
        bldSyncProject();

    MultiFab* vel[MAX_LEV];
    MultiFab* phi[MAX_LEV];
    MultiFab* sig[MAX_LEV];

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
    //
    // Set up outflow bcs.
    //
    if (hasOutFlowBC(phys_bc) && have_divu && do_outflow_bcs)
    {
        set_initial_projection_outflow_bcs(vel,sig,phi,c_lev,cur_divu_time);
    }

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
        }
        //
        // Scale the projection variables.
        //
        scaleVar(sig[lev],1,vel[lev],grids,lev);
    }
    //
    // Setup alias lib.
    //
    const Array<BoxArray>& full_mesh = sync_proj->mesh();
    PArray<MultiFab> u_real[BL_SPACEDIM], p_real(f_lev+1), s_real(f_lev+1);

    int n;
    for (n = 0; n < BL_SPACEDIM; n++) 
        u_real[n].resize(f_lev+1);

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            u_real[n].set(lev, new MultiFab(full_mesh[lev], 1, 1));

            for (MultiFabIterator u_realmfi(u_real[n][lev]); u_realmfi.isValid();
                ++u_realmfi)
            {
                DependentMultiFabIterator velmfi(u_realmfi, *vel[lev]);
                u_realmfi().copy(velmfi(), Xvel+n, 0);
            }
        }
        p_real.set(lev, phi[lev]);
        s_real.set(lev, sig[lev]);
    }
    //
    // Project
    //
    const Real* dx_lev = parent->Geom(f_lev).CellSize();

    if (!have_divu)
    {
        sync_proj->project(u_real, p_real, null_amr_real, s_real,
                           (Real*)dx_lev, proj_tol, c_lev, f_lev,proj_abs_tol);
    } 
    else 
    {
        MultiFab* rhs_cc[MAX_LEV];
        const int nghost = 1; 
        for (lev = c_lev; lev <= f_lev; lev++) 
        {
            AmrLevel&       amr_level = parent->getLevel(lev);
            const BoxArray& grids     = amr_level.boxArray();
            rhs_cc[lev]  = new MultiFab(grids,1,nghost,Fab_allocate);
            MultiFab* rhslev = rhs_cc[lev];
            put_divu_in_cc_rhs(*rhslev,lev,grids,cur_divu_time);
            rhslev->mult(-1.0,0,1,nghost);
            radMult(lev,*rhslev,0); 
        }
        const bool use_u = true;
        PArray<MultiFab> rhs_real(f_lev+1);
        for (lev = c_lev; lev <= f_lev; lev++) 
            rhs_real.set(lev, rhs_cc[lev]);

        sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real,
                                  s_real, use_u, (Real*)dx_lev,
                                  proj_tol, c_lev, f_lev, proj_abs_tol);
        for (lev = c_lev; lev <= f_lev; lev++) 
            delete rhs_cc[lev];
    }
    //
    // Copy and delete u_real, delete s_real if appropriate.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            for (MultiFabIterator u_realmfi(u_real[n][lev]); u_realmfi.isValid();
                ++u_realmfi)
            {
                DependentMultiFabIterator velmfi(u_realmfi, *vel[lev]);
                velmfi().copy(u_realmfi(), 0, Xvel+n);
            }
            delete u_real[n].remove(lev);
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
    //
    // Delete sigma if not used later.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
        delete sig[lev];

    stats.end();
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
    static RunStats stats("sync_project");

    stats.start();

    int lev;
    int f_lev = finest_level;
    if (verbose && ParallelDescriptor::IOProcessor()) 
        cout << "SyncProject: levels = " << c_lev << "  " << f_lev << '\n';
    //
    // Manipulate state + pressure data.
    //
    if (sync_proj == 0)
        bldSyncProject();
    //
    // Gather data.
    //
    MultiFab* vel[MAX_LEV];
    MultiFab* phi[MAX_LEV];
    MultiFab* rhs[MAX_LEV];
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev] = &LevelData[lev].get_new_data(State_Type);
        phi[lev] = &LevelData[lev].get_old_data(Press_Type);
    }

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
            MultiFab& divu_old = amr_level.get_new_data(Divu_Type);
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

            for (MultiFabIterator mfi(*rhslev); mfi.isValid(); ++mfi)
            {
                DependentMultiFabIterator divu_it(mfi,*divu);
                DependentMultiFabIterator dsdt_it(mfi,*dsdt);
                if (!proj_0 && !proj_2) 
                    dsdt_it().minus(divu_it());
                dsdt_it().divide(dt);
                mfi().copy(dsdt_it());
            }

            delete divu;
            delete dsdt;

            rhslev->mult(-1.0,0,1);
            if (CoordSys::IsRZ()) 
                radMult(lev,*rhslev,0);    
        }
    }

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        MultiFab& P_old = LevelData[lev].get_old_data(Press_Type);
        P_old.setVal(0);
    }

    if (hasOutFlowBC(phys_bc) && have_divu && do_outflow_bcs) 
    {
        set_initial_syncproject_outflow_bcs(phi,c_lev,strt_time,dt);
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
    // Convert velocities to accelerations.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        LevelData[lev].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,1);
        LevelData[lev].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,0);
        if (!proj_0 && !proj_2)
        {
            MultiFab& u_o = LevelData[lev].get_old_data(State_Type);
            ConvertUnew(*vel[lev], u_o, dt, LevelData[lev].boxArray());
        }
        else
        {
            const Real dt_inv = 1./dt;
            vel[lev]->mult(dt_inv,0,BL_SPACEDIM,1);
        }
    }
    //
    // Scale initial sync projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        AmrLevel& amr_level = parent->getLevel(lev);
        scaleVar(sig[lev],1,vel[lev],amr_level.boxArray(),lev);
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
        u_real[n].resize(f_lev+1);
    }
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            u_real[n].set(lev, new MultiFab(full_mesh[lev], 1, 1));
            for (MultiFabIterator u_realmfi(u_real[n][lev]); u_realmfi.isValid();
                ++u_realmfi)
            {
                DependentMultiFabIterator velmfi(u_realmfi, *vel[lev]);
                u_realmfi->copy(*velmfi, Xvel+n, 0);
            }
        }
        p_real.set(lev, phi[lev]);
        s_real.set(lev, sig[lev]);
    }

    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        for (lev = f_lev; lev >= c_lev+1; lev--) 
        {
#ifdef BL_USE_HGPROJ_SERIAL
            restrict_level(u_real[n][lev-1],u_real[n][lev],
                           parent->refRatio(lev-1));
#else
            restrict_level(u_real[n][lev-1],u_real[n][lev],parent->refRatio(lev-1),
                           default_restrictor(),level_interface(),0);
#endif
        }
    }

    const Real* dx_lev = parent->Geom(f_lev).CellSize();
    //
    // Project.
    //
    if (!have_divu) 
    {
        //
        // Zero divu only or debugging.
        //
        sync_proj->project(u_real, p_real, null_amr_real, s_real,
                           (Real*)dx_lev, proj_tol, c_lev, f_lev,proj_abs_tol);
    } 
    else 
    {
        //
        // General divu.
        //
        bool use_u = true;
        PArray<MultiFab> rhs_real(f_lev+1);
        for (lev = c_lev; lev <= f_lev; lev++) 
        {
            rhs_real.set(lev, rhs[lev]);
        }
        sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                                  use_u, (Real*)dx_lev,
                                  proj_tol, c_lev, f_lev, proj_abs_tol);
    }
    //
    // Copy and delete u_real, delete s_real if appropriate.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            for (MultiFabIterator u_realmfi(u_real[n][lev]); u_realmfi.isValid();
                ++u_realmfi)
            {
                DependentMultiFabIterator velmfi(u_realmfi, *vel[lev]);
                velmfi().copy(u_realmfi(), 0, Xvel+n);
            }
            delete u_real[n].remove(lev);
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

    stats.end();
}

//
// Compute the node-centered divergence of U.
//

void
Projection::computeDV (MultiFab&       DV,
                       const MultiFab& U,
                       int             src_comp,
                       const Real*     dx,
                       int             is_rz)
{
    const Real      mult    = 1.0;
    const BoxArray& U_boxes = U.boxArray();

    FArrayBox ufab;

    for (MultiFabIterator DVmfi(DV); DVmfi.isValid(); ++DVmfi) 
    {
        DependentMultiFabIterator Umfi(DVmfi, U);

        BL_ASSERT(U_boxes[Umfi.index()] == Umfi.validbox());

        ufab.resize(::grow(Umfi.validbox(),1),BL_SPACEDIM);
        ufab.copy(Umfi(),ufab.box(),src_comp,ufab.box(),0,BL_SPACEDIM);

        Box ndbox = ::surroundingNodes(Umfi.validbox());

        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();
        const int* ulo  = ufab.box().loVect();
        const int* uhi  = ufab.box().hiVect();

        FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
                    DVmfi().dataPtr(),ARLIM(ndlo),ARLIM(ndhi),
                    ndlo,ndhi,dx,&mult,&is_rz);
    }
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
    int             lowfix = (isrz==1 && bcxlo!=FOEXTRAP && bcxlo!=HOEXTRAP);
    int             hifix  = (isrz==1 && bcxhi!=FOEXTRAP && bcxhi!=HOEXTRAP);
    const Real*     dx     = geom.CellSize();
    Real            hx     = dx[0];
    const Box&      domain = geom.Domain();
    const int*      domlo  = domain.loVect();
    const int*      domhi  = domain.hiVect();

    NavierStokes* ns = dynamic_cast<NavierStokes*>(&parent->getLevel(level));

    BL_ASSERT(!(ns == 0));

    MultiFab* divu = ns->getDivCond(1,time);

    for (MultiFabIterator rhsmfi(rhs); rhsmfi.isValid(); ++rhsmfi)
    {
        DependentMultiFabIterator divumfi(rhsmfi,*divu);

        DEF_CLIMITS(divumfi(),divudat,divulo,divuhi);
        DEF_LIMITS(rhsmfi(),rhsdat,rhslo,rhshi);
        Array<Real> rcen(divumfi().box().length(0),1.0);
        if (isrz == 1) 
            geom.GetCellLoc(rcen,divumfi().box(),0);

        FORT_HGC2N(&nghost,ARLIM(divulo),ARLIM(divuhi),divudat,
                   rcen.dataPtr(), ARLIM(rhslo),ARLIM(rhshi),rhsdat,
                   domlo,domhi,lowfix,hifix,&hx,&isrz);
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

    for (MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi, *divu);

        mfi().copy(dmfi());
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

    geom.FillPeriodicBoundary(psi,0,nvar,false,false);
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
    for (MultiFabIterator Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        DependentMultiFabIterator Unewmfi(Uoldmfi, Unew);

        BL_ASSERT(grids[Uoldmfi.index()] == Uoldmfi.validbox());

        UnConvertUnew( Uoldmfi(), alpha, Unewmfi(), Uoldmfi.validbox() );
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
    for (MultiFabIterator Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        DependentMultiFabIterator Unewmfi(Uoldmfi, Unew);

        BL_ASSERT(grids[Uoldmfi.index()] == Uoldmfi.validbox());

        ConvertUnew( Unewmfi(), Uoldmfi(), alpha, Uoldmfi.validbox() );
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
    for (MultiFabIterator Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        DependentMultiFabIterator Unewmfi(Uoldmfi, Unew);

        BL_ASSERT(grids[Uoldmfi.index()] == Uoldmfi.validbox());

        UpdateArg1(Unewmfi(),alpha,Uoldmfi(),nvar,Uoldmfi.validbox(),ngrow);
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

    Box        b  = ::grow(grd,ngrow);
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
    for (MultiFabIterator pmfi(p); pmfi.isValid(); ++pmfi) 
    {
        DependentMultiFabIterator phimfi(pmfi, phi);
        pmfi().plus(phimfi());
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

    for (MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi)
    {
        DependentMultiFabIterator P_oldmfi(P_newmfi, P_old);

        UpdateArg1(P_newmfi(),1.0/dt,P_oldmfi(),1,grids[P_newmfi.index()],1);

        P_oldmfi().setVal(BogusValue);
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

    int n[BL_SPACEDIM];

    for (int i = 0; i < BL_SPACEDIM; i++)
        n[i] = ngrow;

    for (MultiFabIterator mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
    {
        BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

        const int* lo = mfmfi.validbox().loVect();
        const int* hi = mfmfi.validbox().hiVect();
        Real* dat     = mfmfi().dataPtr(comp);
        Real* rad     = &radius[level][mfmfi.index()];

        FORT_RADMPY(dat,ARLIM(lo),ARLIM(hi),&ngrow,rad,&nr,n);
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

    int n[BL_SPACEDIM];

    for (int i = 0; i < BL_SPACEDIM; i++)
        n[i] = ngrow;

    for (MultiFabIterator mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
    {
        BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

        const int* lo  = mfmfi.validbox().loVect();
        const int* hi  = mfmfi.validbox().hiVect();
        Real*      dat = mfmfi().dataPtr(comp);
        Real*      rad = &radius[level][mfmfi.index()];

        FORT_RADDIV(dat,ARLIM(lo),ARLIM(hi),&ngrow,rad,&nr,n);
    }
}

void
Projection::set_level_projector_outflow_bcs (int       level,
                                             Real      mid_time,
                                             MultiFab& phi,
                                             MultiFab& vel,
                                             MultiFab& divu)
{
    //
    // Do as set_initial_projection_outflow_bcs, except this one for single
    // specified level.  In this case, since the finest level has not been
    // advanced to the new time, the data there is not yet available.  We
    // contruct the dirichlet phi values for the outflow boundary at the
    // requested level, fill-patching in coarser data, which presumably has
    // been advanced up to or past new_time.
    //
    // FIXME: the three outflow boundary routines should be collapsed into
    //        one, since they all do roughly the same thing.
    //
#if (BL_SPACEDIM == 2)
    //
    // Make sure out flow only occurs at yhi faces.
    //
    const Orientation outFace(1, Orientation::high);
    bool hasOutFlow;
    Orientation _outFace;
    getOutFlowFace(hasOutFlow,_outFace,phys_bc);
    BL_ASSERT(_outFace == outFace);

    const Real* dx         = parent->Geom(level).CellSize();
    const Box& domain      = parent->Geom(level).Domain();
    const int outDir       = 1;
    const int ccStripWidth = 3;
    const int ncStripWidth = 1;
    const Box state_strip  = Box(::adjCellHi(domain,outDir,ccStripWidth)).shift(outDir,-ccStripWidth+1);
    const Box phi_strip    = ::surroundingNodes(bdryHi(domain,outDir,ncStripWidth));
    //
    // Only do this if the entire outflow face is covered by boxes at this level
    // (skip if doesnt touch, and bomb if only partially covered)
    //
    const BoxArray& grids = parent->getLevel(level).boxArray();
    const Box valid_state_strip = state_strip & domain;
    const BoxArray uncovered_outflow_ba = ::complementIn(valid_state_strip,grids);

    BL_ASSERT( !(uncovered_outflow_ba.ready() &&
              ::intersect(grids,valid_state_strip).ready()) );

    if ( !(uncovered_outflow_ba.ready()) )
    {
        const int nGrow        = 0;
        const int nCompPhi     = 1;
        const int srcCompRho = Density, nCompRho = 1;
        const BoxArray state_strip_ba(&state_strip,1);
        const BoxArray phi_strip_ba(&phi_strip,1);

        MultiFab cc_divu(state_strip_ba, 1, nGrow, Fab_allocate);
        MultiFab cc_vel(state_strip_ba, BL_SPACEDIM, nGrow, Fab_allocate);
        MultiFab phi_fine_strip(phi_strip_ba, nCompPhi, nGrow, Fab_allocate);
        //
        // Fill Fab-like MultiFab of divu data, and initialize phi
        //
        cc_divu.copy(divu);
        cc_vel.copy(vel);
        phi_fine_strip.setVal(0);
        //
        // Make r_i needed in HGPHIBC (set = 1 if cartesian).
        //
        Box region = Box(::adjCellHi(domain,outDir,1)).shift(outDir, -1);
        Array<Real> rcen(region.length(0), 1.0);
        if (CoordSys::IsRZ() == 1) 
            parent->Geom(level).GetCellLoc(rcen, region, 0);
    
        FillPatchIterator rhoFpi(LevelData[level],cc_divu,nGrow,
                                 mid_time,State_Type,srcCompRho,nCompRho);

        for ( ; rhoFpi.isValid(); ++rhoFpi)
        {
            DependentMultiFabIterator phimfi(rhoFpi, phi_fine_strip);
            DependentMultiFabIterator vel_mfi(rhoFpi,cc_vel);
            DependentMultiFabIterator divu_mfi(rhoFpi,cc_divu);
            //
            // Fill phi_fine_strip with boundary cell values for phi, then copy
            // into arg data Note: Though this looks like a distributed
            // operation, the MultiFab is built on a single box.  This is
            // necessary currently, since FORT_HGPHIBC requires a slice across
            // the entire domain.
            //
            const int isPeriodicInX = parent->Geom(level).isPeriodic(0);

            FORT_HGPHIBC(ARLIM(vel_mfi().loVect()), ARLIM(vel_mfi().hiVect()), vel_mfi().dataPtr(),
                         ARLIM(divu_mfi().loVect()),ARLIM(divu_mfi().hiVect()),divu_mfi().dataPtr(),
                         ARLIM(rhoFpi().loVect()),  ARLIM(rhoFpi().hiVect()),  rhoFpi().dataPtr(),
                         ARLIM(region.loVect()),    ARLIM(region.hiVect()),    rcen.dataPtr(),   &dx[0],
                         ARLIM(phimfi().loVect()),  ARLIM(phimfi().hiVect()),  phimfi().dataPtr(),
                         &isPeriodicInX);
        }

        phi.copy(phi_fine_strip);
    }
#else
    // check to see if divu == 0 near outflow.  If it isn't, then abort.
    REAL divu_norm = checkDivU(parent,LevelData,phys_bc,
			       level,level,mid_time);
    if (divu_norm > 1.0e-7) {
      cout << "divu_norm = " << divu_norm << endl;
      BoxLib::Error("outflow bc for divu != 0 not implemented in 3D");
    }
#endif
}

void Projection::set_initial_projection_outflow_bcs (MultiFab** vel,
                                                     MultiFab** sig,
                                                     MultiFab** phi,
                                                     int        c_lev,
                                                     Real       cur_divu_time)
{
#if (BL_SPACEDIM == 2)
    //
    // Warning: This code looks about right, but hasn't really been tested.
    //
    // what we do
    //  at the finest level covering the entire outflow boundary
    //  1) we define 3 cell-wide strips for
    //     rho, S, U at the top of the domain
    //     The strips are three cells wide for convenience
    //  2) we define a 1 node-wide for phi at the top of the domain
    //  3) we fill the rho, S, and U strips
    //  4) we compute phi on the phi strip with HGPHIBC
    //  5) we copy the phi strip into the phi multifab
    //  we then "putdown" onto coarser levels
    //
    // FIXME: vel and sig unused here.  We needed a filpatch operation, so pulled
    //        data from the state (presumably is where vel, sig from anyway).  sig
    //        is poorly named, and we need other stuff as well, so we should
    //        probably just remove sig and vel from args
    //
    // make sure out flow only occurs at yhi faces
    const Orientation outFace(1, Orientation::high);
    bool hasOutFlow;
    Orientation _outFace;
    getOutFlowFace(hasOutFlow,_outFace,phys_bc);
    BL_ASSERT(_outFace == outFace);

    BL_ASSERT(c_lev == 0);
    //
    // Get 3-wide cc box, state_strip, along top, incl. 2 rows of int and
    // 1 row of ghosts.
    // Get 1-wide nc box, phi_strip, along top, grown out by one
    // tangential to face.
    //
    const int outDir       = 1;
    const int ccStripWidth = 3;
    const int ncStripWidth = 1;
    //
    // Determine the finest level such that the entire outflow face is covered
    // by boxes at this level (skip if doesnt touch, and bomb if only partially
    // covered)
    //
    Box state_strip, domain;
    int f_lev = finest_level;
    for ( ; f_lev>=c_lev; --f_lev)
    {
        domain = parent->Geom(f_lev).Domain();
        state_strip
            = Box(::adjCellHi(domain,outDir,ccStripWidth)).shift(outDir,-ccStripWidth+1);

        const BoxArray& Lgrids = parent->getLevel(f_lev).boxArray();
        const Box valid_state_strip = state_strip & domain;
        const BoxArray uncovered_outflow_ba = ::complementIn(valid_state_strip,Lgrids);
        BL_ASSERT( !(uncovered_outflow_ba.ready() &&
                  ::intersect(Lgrids,valid_state_strip).ready()) );
        if ( !(uncovered_outflow_ba.ready()) )
            break;
    }
    const Real* dx = parent->Geom(f_lev).CellSize();
    const Box phi_strip = ::surroundingNodes(bdryHi(domain,outDir,ncStripWidth));

    const int   nGrow    = 0;
    const int   nCompPhi = 1;
    const int srcCompRho = Density, nCompRho = 1;
    const int srcCompVel = Xvel,    nCompVel = BL_SPACEDIM;
    const int srcCompDivu = 0,      nCompDivu = 1;
    
    BoxArray state_strip_ba(&state_strip,1);
    BoxArray phi_strip_ba(&phi_strip,1);
    
    MultiFab cc_MultiFab(state_strip_ba, 1, nGrow, Fab_noallocate);
    MultiFab phi_fine_strip(phi_strip_ba, nCompPhi, nGrow, Fab_allocate);
    
    phi_fine_strip.setVal(0);
    //
    // Make r_i needed in HGPHIBC (set = 1 if cartesian).
    //
    Box region = Box(::adjCellHi(domain,outDir,1)).shift(outDir, -1);
    Array<Real> rcen(region.length(0), 1.0);
    if (CoordSys::IsRZ() == 1) 
        parent->Geom(f_lev).GetCellLoc(rcen, region, 0);
    
    int Divu_Type, Divu;
    FArrayBox rho;
    if (!LevelData[c_lev].isStateVariable("divu", Divu_Type, Divu)) 
        BoxLib::Error("Projection::set_initial_projection_outflow_bcs(): no Divu");
    
    FillPatchIterator rhoFpi(LevelData[f_lev],cc_MultiFab);
    if (rho_wgt_vel_proj)
        rhoFpi.Initialize(nGrow,cur_divu_time,State_Type,srcCompRho,nCompRho);
    
    FillPatchIterator velFpi(LevelData[f_lev],cc_MultiFab,nGrow,
                             cur_divu_time,State_Type,srcCompVel,nCompVel);
    
    FillPatchIterator divuFpi(LevelData[f_lev],cc_MultiFab,nGrow,
                              cur_divu_time,Divu_Type,srcCompDivu,nCompDivu);

    for ( ;velFpi.isValid() && divuFpi.isValid(); ++velFpi, ++divuFpi)
    {
        DependentMultiFabIterator phimfi(rhoFpi, phi_fine_strip);
        rho.resize(velFpi.validbox(), nCompRho);
        if (rho_wgt_vel_proj && rhoFpi.isValid())
            rho.copy(rhoFpi());
        else
            rho.setVal(1.0);
        //
        // Fill phi_fine_strip with boundary cell values for phi, then
        // copy into arg data.  Note: Though this looks like a distributed
        // operation, the MultiFab is built on a single box...this is
        // necesary currently, since FORT_HGPHIBC requires a slice across
        // the entire domain.
        //
        const int isPeriodicInX = parent->Geom(f_lev).isPeriodic(0);
        
        FORT_HGPHIBC(ARLIM(velFpi().loVect()),  ARLIM(velFpi().hiVect()),  velFpi().dataPtr(),
                     ARLIM(divuFpi().loVect()), ARLIM(divuFpi().hiVect()), divuFpi().dataPtr(),
                     ARLIM(rho.loVect()),       ARLIM(rho.hiVect()),       rho.dataPtr(),
                     ARLIM(region.loVect()),    ARLIM(region.hiVect()),    rcen.dataPtr(),    &dx[0],
                     ARLIM(phimfi().loVect()),  ARLIM(phimfi().hiVect()),  phimfi().dataPtr(),
                     &isPeriodicInX);

        if (rho_wgt_vel_proj)
            ++rhoFpi;
    }
    
    phi[f_lev]->copy(phi_fine_strip);
    
    IntVect ratio = IntVect::TheUnitVector();
    
    for (int lev = f_lev-1; lev >= c_lev; lev--) 
    {
        ratio *= parent->refRatio(lev);
        const Box& domainC = parent->Geom(lev).Domain();
        Box top_phiC_strip = ::surroundingNodes(bdryHi(domainC,outDir,ncStripWidth));
        BoxArray top_phiC_strip_ba(&top_phiC_strip,1);
        MultiFab phi_crse_strip(top_phiC_strip_ba,nCompPhi,nGrow);
        phi_crse_strip.setVal(0);
        
        for (MultiFabIterator finemfi(phi_fine_strip); finemfi.isValid(); ++finemfi)
        {
            DependentMultiFabIterator crsemfi(finemfi, phi_crse_strip);
            Box ovlp = Box(finemfi.validbox()).coarsen(ratio) & crsemfi.validbox();
            FORT_PUTDOWN (crsemfi().dataPtr(),ARLIM(crsemfi().loVect()),ARLIM(crsemfi().hiVect()),
                          finemfi().dataPtr(),ARLIM(finemfi().loVect()),ARLIM(finemfi().hiVect()),
                          ovlp.loVect(),ovlp.hiVect(),ratio.getVect());
        }
        phi[lev]->copy(phi_crse_strip);
    }
#else
    // check to see if divu == 0 near outflow.  If it isn't, then abort.
    REAL divu_norm = checkDivU(parent,LevelData,phys_bc,finest_level,c_lev,cur_divu_time);
    if (divu_norm > 1.0e-7) {
      cout << "divu_norm = " << divu_norm << endl;
      BoxLib::Error("outflow bc for divu != 0 not implemented in 3D");
    }
#endif
}

void
Projection::set_initial_syncproject_outflow_bcs (MultiFab** phi, 
                                                 int        c_lev,
                                                 Real       start_time,
                                                 Real       dt)
{
#if (BL_SPACEDIM == 2)
    //
    // Warning: This code looks about right, but hasn't really been tested.
    //
    // What we do is similar to what we do in set_initial_projection_outflow_bcs
    // except that we work with time derivatives of U and S.
    //

    // make sure out flow only occurs at yhi faces
    const Orientation outFace(1, Orientation::high);
    bool hasOutFlow;
    Orientation _outFace;
    getOutFlowFace(hasOutFlow,_outFace,phys_bc);
    BL_ASSERT(_outFace == outFace);

    const int f_lev = finest_level;

    BL_ASSERT(c_lev == 0);

    const Real* dx           = parent->Geom(f_lev).CellSize();
    const Box&  domain       = parent->Geom(f_lev).Domain();
    const int   outDir       = 1;
    const int   ccStripWidth = 3;
    const int   ncStripWidth = 1;
    Box state_strip =
        ::adjCellHi(domain,outDir,ccStripWidth).shift(outDir,-ccStripWidth+1);
    Box phi_strip = 
        ::surroundingNodes(bdryHi(domain,outDir,ncStripWidth));
    const int   nGrow        = 0;
    const int   nCompPhi     = 1;

    const int srcCompRho = Density, nCompRho = 1;
    const int srcCompVel = Xvel,    nCompVel = BL_SPACEDIM;
    const int srcCompDivu = 0,      nCompDivu = 1;

    BoxArray state_strip_ba(&state_strip,1);
    BoxArray phi_strip_ba(&phi_strip,1);

    MultiFab cc_MultiFab(state_strip_ba, 1, nGrow, Fab_noallocate);
    MultiFab phi_fine_strip(phi_strip_ba, nCompPhi, nGrow, Fab_allocate);

    phi_fine_strip.setVal(0);
    //
    // Make r_i needed in HGPHIBC (set = 1 if cartesian).
    //
    Box region = Box(::adjCellHi(domain,outDir,1)).shift(outDir, -1);
    Array<Real> rcen(region.length(0), 1.0);
    if (CoordSys::IsRZ() == 1) 
        parent->Geom(f_lev).GetCellLoc(rcen, region, 0);
    
    int Divu_Type, Divu;
    if (!LevelData[c_lev].isStateVariable("divu", Divu_Type, Divu)) 
        BoxLib::Error("Projection::set_syncproject_outflow_bcs(): Divu not found");

    FArrayBox rhonph, dudt, dsdt;
    
    FillPatchIterator rhoOldFpi(LevelData[f_lev],cc_MultiFab,nGrow,
                                start_time,State_Type,srcCompRho,nCompRho);
	
    FillPatchIterator rhoNewFpi(LevelData[f_lev],cc_MultiFab,nGrow,
                                start_time+dt,State_Type,srcCompRho,nCompRho);
	
    FillPatchIterator velOldFpi(LevelData[f_lev],cc_MultiFab,nGrow,
                                start_time,State_Type,srcCompVel,nCompVel);
	
    FillPatchIterator velNewFpi(LevelData[f_lev],cc_MultiFab,nGrow,
                                start_time+dt,State_Type,srcCompVel,nCompVel);
	
    FillPatchIterator divuOldFpi(LevelData[f_lev],cc_MultiFab,nGrow,
                                 start_time,Divu_Type,srcCompDivu,nCompDivu);
	
    FillPatchIterator divuNewFpi(LevelData[f_lev],cc_MultiFab,nGrow,
                                 start_time+dt,Divu_Type,srcCompDivu,nCompDivu);
	
    for ( ;rhoOldFpi.isValid()  &&  velOldFpi.isValid()  &&  divuOldFpi.isValid() &&
              rhoNewFpi.isValid()  &&  velNewFpi.isValid()  &&  divuNewFpi.isValid();
          ++rhoOldFpi, ++velOldFpi, ++divuOldFpi, ++rhoNewFpi, ++velNewFpi, ++divuNewFpi)
    {
        DependentMultiFabIterator phimfi(rhoOldFpi, phi_fine_strip);
        //
        // Make rhonph, du/dt, and dsdt.
        //
        BL_ASSERT(rhoOldFpi.validbox() == rhoNewFpi.validbox());
        rhonph.resize(rhoOldFpi.validbox(),nCompRho);
        rhonph.copy(rhoNewFpi());
        rhonph.plus(rhoOldFpi());
        rhonph.divide(2.0);

        BL_ASSERT(divuOldFpi.validbox() == divuNewFpi.validbox());
        dudt.resize(velOldFpi.validbox(),nCompVel);
        dudt.copy(velNewFpi());
        if (!proj_0 && !proj_2)
            dudt.minus(velOldFpi());
        dudt.divide(dt);
	    
        BL_ASSERT(velOldFpi.validbox() == velNewFpi.validbox());
        dsdt.resize(divuOldFpi.validbox(),nCompDivu);
        dsdt.copy(divuNewFpi());
        if (!proj_0 && !proj_2)
            dsdt.minus(divuOldFpi());
        dsdt.divide(dt);
        //
        // Fill phi_fine_strip with boundary cell values for phi, then
        // copy into arg data.  Note: Though this looks like a distributed
        // operation, the MultiFab is built on a single box...this is
        // necesary currently, since FORT_HGPHIBC requires a slice across
        // the entire domain
        //
        const int isPeriodicInX = parent->Geom(f_lev).isPeriodic(0);

        FORT_HGPHIBC(ARLIM(dudt.loVect()),     ARLIM(dudt.hiVect()),     dudt.dataPtr(),
                     ARLIM(dsdt.loVect()),     ARLIM(dsdt.hiVect()),     dsdt.dataPtr(),
                     ARLIM(rhonph.loVect()),   ARLIM(rhonph.hiVect()),   rhonph.dataPtr(),
                     ARLIM(region.loVect()),   ARLIM(region.hiVect()),   rcen.dataPtr(),     &dx[0],
                     ARLIM(phimfi().loVect()), ARLIM(phimfi().hiVect()), phimfi().dataPtr(),
                     &isPeriodicInX);
    }

    phi[f_lev]->copy(phi_fine_strip);
    //
    // Put down to coarser levels.
    //
    IntVect ratio = IntVect::TheUnitVector();
    for (int lev = f_lev-1; lev >= c_lev; lev--) 
    {
        ratio *= parent->refRatio(lev);
        const Box& domainC = parent->Geom(lev).Domain();
        Box top_phiC_strip = ::surroundingNodes(bdryHi(domainC,outDir,ncStripWidth));
	    
        BoxArray top_phiC_strip_ba(&top_phiC_strip,1);
        MultiFab phi_crse_strip(top_phiC_strip_ba, nCompPhi, nGrow);
        phi_crse_strip.setVal(0);
	    
        for (MultiFabIterator finemfi(phi_fine_strip); finemfi.isValid(); ++finemfi)
        {
            DependentMultiFabIterator crsemfi(finemfi, phi_crse_strip);
            Box ovlp = Box(finemfi.validbox()).coarsen(ratio) & crsemfi.validbox();
            FORT_PUTDOWN (crsemfi().dataPtr(),ARLIM(crsemfi().loVect()),ARLIM(crsemfi().hiVect()),
                          finemfi().dataPtr(),ARLIM(finemfi().loVect()),ARLIM(finemfi().hiVect()),
                          ovlp.loVect(),ovlp.hiVect(),ratio.getVect());
        }
        phi[lev]->copy(phi_crse_strip);
    }
#else
    // check to see if divu == 0 near outflow.  If it isn't, then abort.
    REAL oldDivu_norm = checkDivU(parent,LevelData,phys_bc,
				  finest_level,c_lev,start_time);
    // check on both norms?
    REAL newDivu_norm = checkDivU(parent,LevelData,phys_bc,
				  finest_level,c_lev,start_time+dt);
    if (oldDivu_norm > 1.0e-7) {
      cout << "divu_norm = " << oldDivu_norm << endl;
      BoxLib::Error("outflow bc for divu != 0 not implemented in 3D");
    }
#endif
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
        cout << "initialVorticityProject: levels = "
             << c_lev
             << "  "
             << f_lev << endl;
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

    projector_bndry = new inviscid_fluid_boundary_class(proj_bc);

    bldSyncProject();

    MultiFab* vel[MAX_LEV];
    MultiFab* phi[MAX_LEV];
    MultiFab* sig[MAX_LEV];

    for (int lev = c_lev; lev <= f_lev; lev++)
    {
        MultiFab& P_new = LevelData[lev].get_new_data(Press_Type);

        const int nghost = 1;

        sig[lev] = new MultiFab(LevelData[lev].boxArray(),1,nghost);
        sig[lev]->setVal(1,nghost);
        phi[lev] = new MultiFab(P_new.boxArray(),1,nghost);
        phi[lev]->setVal(0,nghost);
    }
    //
    // Set up outflow bcs.
    //
    const Array<BoxArray>& full_mesh = sync_proj->mesh();
    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> p_real(f_lev+1), s_real(f_lev+1);
    PArray<MultiFab> rhs_real(f_lev+1);

    for (int n = 0; n < BL_SPACEDIM; n++)
        u_real[n].resize(f_lev+1);

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

        for (MultiFabIterator mfi(rhs_real[lev]); mfi.isValid(); ++mfi)
        {
            DependentMultiFabIterator dmfi(mfi,P_new);

            mfi().setVal(0);
            mfi().copy(dmfi(), 0, 0);
        }

        p_real.set(lev, phi[lev]);
        s_real.set(lev, sig[lev]);
    }
    //
    // Project.
    //
    const Real* dx_lev = parent->Geom(f_lev).CellSize();
    const int   use_u  = 0;
    sync_proj->manual_project(u_real,p_real,rhs_real,null_amr_real,s_real,
                              use_u,(Real*)dx_lev,
                              proj_tol,c_lev,f_lev,proj_abs_tol);
    //
    // Copy and delete u_real, delete s_real if appropriate.
    //
    for (int lev = c_lev; lev <= f_lev; lev++)
    {
        vel[lev] = &LevelData[lev].get_new_data(State_Type);

        int n    = 0;
        int ninv = 1;

        for (MultiFabIterator mfi(*vel[lev]); mfi.isValid(); ++mfi)
        {
            DependentMultiFabIterator dmfi(mfi,u_real[n][lev]);

            mfi().copy(dmfi(), 0, Xvel+ninv);
            mfi().mult(-1, Xvel+ninv, 1);
        }

        delete u_real[n].remove(lev);

        n    = 1;
        ninv = 0;

        for (MultiFabIterator mfi(*vel[lev]); mfi.isValid(); ++mfi)
        {
            DependentMultiFabIterator dmfi(mfi,u_real[n][lev]);

            mfi().copy(dmfi(), 0, Xvel+ninv);
        }

        delete u_real[n].remove(lev);

        delete rhs_real.remove(lev);
    }

    for (int lev = c_lev; lev <= f_lev; lev++)
        delete sig[lev];
    //
    // Reset the boundary conditions for all the other projections.
    //
    delete sync_proj;

    sync_proj = 0;
#else
    BoxLib::Error("Projection::initialVorticityProject(): not implented yet for 3D");
#endif
}

