//BL_COPYRIGHT_NOTICE

//
// $Id: MacProj.cpp,v 1.42 1999-02-26 18:53:12 propp Exp $
//

#include <Misc.H>
#include <LO_BCTYPES.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <MacProj.H>
#include <MacBndry.H>
#include <MacOpMacDrivers.H>
#include <NavierStokes.H>
#include <MACPROJ_F.H>

const char NL = '\n';

#ifndef _NavierStokes_H_
enum StateType {State_Type=0, Press_Type};
#  if (BL_SPACEDIM == 2)
enum StateNames  { Xvel=0, Yvel, Density};
#  else
enum StateNames  { Xvel=0, Yvel, Zvel, Density};
#  endif
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define GEOM_GROW 1
#define HYP_GROW 3

int  MacProj::verbose          = 0;
int  MacProj::use_cg_solve     = 0;
Real MacProj::mac_tol          = 1.0e-12;
Real MacProj::mac_abs_tol      = 1.0e-16;
Real MacProj::mac_sync_tol     = 1.0e-8;
int  MacProj::do_outflow_bcs   = 1;
int  MacProj::fix_mac_sync_rhs = 0;

//
// Setup functions follow
//

MacProj::MacProj (Amr*   _parent,
                  int    _finest_level,
                  BCRec* _phys_bc,
                  int    _radius_grow)
  :
    parent(_parent),
    finest_level(_finest_level),
    phys_bc(_phys_bc), 
    radius_grow(_radius_grow), 
    LevelData(_finest_level+1),
    phi_bcs(_finest_level+1),
    mac_phi_crse(_finest_level+1, PArrayManage),
    mac_reg(_finest_level+1, PArrayManage),
    volume(_finest_level+1), area(_finest_level+1),
    radius(_finest_level+1) 
{
    read_params();

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "Creating mac_projector\n";
    }

    finest_level_allocated = finest_level;
}

MacProj::~MacProj () {}

void
MacProj::read_params ()
{
    //
    // Read parameters from input file and command line.
    //
    ParmParse pp("mac");

    pp.query( "v",                verbose          );

    pp.query( "mac_tol",          mac_tol          );
    pp.query( "mac_sync_tol",     mac_sync_tol     );
    pp.query( "use_cg_solve",     use_cg_solve     );
    pp.query( "mac_abs_tol",      mac_abs_tol      );
    pp.query( "do_outflow_bcs",   do_outflow_bcs   );
    pp.query( "fix_mac_sync_rhs", fix_mac_sync_rhs );
}

void
MacProj::install_level (int           level,
                        AmrLevel*     level_data,
                        MultiFab&     _volume,
                        MultiFab*     _area,
                        PArray<Real>* _radius )
{
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "Installing MacProj level " << level << NL;
    }

    if (parent->finestLevel() < finest_level)
    {
        for (int lev = parent->finestLevel() + 1; lev <= finest_level; lev++)
        {
            mac_reg.clear(lev);
        }
    }

    finest_level = parent->finestLevel();

    if (level > finest_level_allocated)
    {
        finest_level_allocated = finest_level;
        LevelData.resize(finest_level+1);
        phi_bcs.resize(finest_level+1);
        mac_phi_crse.resize(finest_level+1);
        mac_reg.resize(finest_level+1);
        volume.resize(finest_level+1);
        area.resize(finest_level+1);
        radius.resize(finest_level+1);
    }

    LevelData.clear(level);
    LevelData.set(level, level_data);
    volume.clear(level);
    volume.set(level, &_volume);
    area.set(level, _area);
    radius.clear(level);
    radius.set(level, _radius);

    BuildPhiBC(level);

    if (level > 0)
    {
        const BoxArray& grids = LevelData[level].boxArray();
        mac_reg.clear(level);
        mac_reg.set(level, new FluxRegister(grids,parent->refRatio(level-1),
                                            level,1));
    }
}

void
MacProj::BuildPhiBC (int level)
{
    const BoxArray& grids = LevelData[level].boxArray();
    const Geometry& geom  = parent->Geom(level);
    const int ngrds       = grids.length();
    phi_bcs[level].resize(ngrds);
    const Box& domain     = geom.Domain();
    const int* domlo      = domain.loVect();
    const int* domhi      = domain.hiVect();
    const int* phys_lo    = phys_bc->lo();
    const int* phys_hi    = phys_bc->hi();

    for (int i = 0; i < ngrds; i++)
    {
        BCRec& bc     = phi_bcs[level][i];
        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (lo[dir] == domlo[dir])
            {
                if (phys_lo[dir] == Outflow)
                    bc.setLo(dir,LO_DIRICHLET);
                else
                    bc.setLo(dir,LO_NEUMANN);
            }
            else
            {
                bc.setLo(dir,LO_DIRICHLET);
            }
            if (hi[dir] == domhi[dir])
            {
                if (phys_hi[dir] == Outflow)
                    bc.setHi(dir,LO_DIRICHLET);
                else
                    bc.setHi(dir,LO_NEUMANN);
            }
            else
            {
                bc.setHi(dir,LO_DIRICHLET);
            }
        }
    }
}

void
MacProj::setup (int level)
{
    if (level < parent->maxLevel())
    {
        if (!mac_phi_crse.defined(level))
        {
            const BoxArray& grids = LevelData[level].boxArray();
            mac_phi_crse.set(level, new MultiFab(grids,1,1,Fab_allocate));
            mac_phi_crse[level].setVal(0.0);
        }
    }
}

void
MacProj::cleanup (int level)
{
    if (level < parent->maxLevel())
    {
        mac_phi_crse.clear(level);
    }
}

//
// Projection functions follow ...
//

static
bool
grids_on_side_of_domain (const BoxArray& grids,
                        const Box&      domain,
			const Orientation& outFace)
{
  const int idir = outFace.coordDir();

  if (outFace.isLow()) {
    for (int igrid = 0; igrid < grids.length(); igrid++) { 
      if (grids[igrid].smallEnd(idir) == domain.smallEnd(idir)) { 
	return true;
      }
    }
  }
  
  if (outFace.isHigh()) {
    for (int igrid = 0; igrid < grids.length(); igrid++) {
      if (grids[igrid].bigEnd(idir) == domain.bigEnd(idir)) {
	return true;
      }
    }
  }

  return false;
}

//
// Compute the level advance mac projection.
//

void
MacProj::mac_project (int             level,
                      MultiFab*       u_mac,
                      MultiFab&       S,
                      Real            dt,
                      Real            time,
                      const MultiFab& divu,
                      int             have_divu)
{
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "... mac_project at level " << level << NL;
    }
    const BoxArray& grids = LevelData[level].boxArray();
    const Geometry& geom  = parent->Geom(level);
    const Real* dx        = geom.CellSize();
    const int max_level   = parent->maxLevel();
    MultiFab* mac_phi     = 0;
    IntVect crse_ratio    = (level > 0) ? 
        parent->refRatio(level-1) : IntVect::TheZeroVector();
    //
    // If finest level possible no need to make permanent mac_phi for bcs.
    //
    if (level == max_level)
    {
        mac_phi = new MultiFab(grids,1,1,Fab_allocate);
    }
    else
    {
        mac_phi = &mac_phi_crse[level];
    }

    mac_phi->setVal(0.0);

#if (BL_SPACEDIM == 2)
    int outflow_at_top = phys_bc->lo(0) != Outflow && phys_bc->lo(1) != Outflow && 
        phys_bc->hi(0) != Outflow && phys_bc->hi(1) == Outflow;
    const Orientation outFace(1,Orientation::high);
    if (outflow_at_top &&
        have_divu      &&
        do_outflow_bcs &&
        grids_on_side_of_domain(grids,geom.Domain(),outFace))
    {
        set_outflow_bcs(level, mac_phi, u_mac, S, divu);
    }
#endif

    //
    // Store the Dirichlet boundary condition for mac_phi in mac_bndry.
    //
    MacBndry mac_bndry(grids,1,geom);
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;
    if (level == 0)
    {
        mac_bndry.setBndryValues(*mac_phi,src_comp,dest_comp,num_comp,*phys_bc);
    }
    else
    {
        MultiFab& CPhi = mac_phi_crse[level-1];
        BoxArray crse_boxes(grids);
        crse_boxes.coarsen(crse_ratio);
        const int in_rad     = 0;
        const int out_rad    = 1;
        const int extent_rad = 1;
        BndryRegister crse_br(crse_boxes,in_rad,out_rad,extent_rad,num_comp);
        crse_br.copyFrom(CPhi,extent_rad,src_comp,dest_comp,num_comp);

        mac_bndry.setBndryValues(crse_br,src_comp,*mac_phi,src_comp,
                                 dest_comp,num_comp,crse_ratio,*phys_bc);
    }
    //
    // Compute the nondivergent velocities, by creating the linop
    // and multigrid operator appropriate for the solved system.
    //
    // Initialize the rhs with divu.
    //
    const Real rhs_scale = 2.0/dt;
    MultiFab Rhs(grids,1,0,Fab_allocate);
    Rhs.copy(divu);

    mac_level_driver(mac_bndry, grids, use_cg_solve, level, Density,
                     dx, dt, mac_tol, mac_abs_tol, rhs_scale, 
                     area[level], volume[level], S, Rhs, u_mac, mac_phi);
    //
    // Test that u_mac is divergence free
    //
    if (verbose)
        check_div_cond(level, u_mac);
    //
    // Store advection velocities in mac registers at crse/fine boundaries.
    //
    // Initialize advection velocity registers with coarse grid velocity.
    //
    if (level < finest_level)
    {
        FluxRegister& mr = mac_reg[level+1];
        mr.setVal(0.0);
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            mr.CrseInit(u_mac[dir],area[level][dir],dir,0,0,1,-1.0);
        }

        Real sumreg =  mr.SumReg(0);
        if (verbose && ParallelDescriptor::IOProcessor())
        {
            cout << "LEVEL "
                 << level
                 << " MACREG: CrseInit sum = "
                 << sumreg << endl;
        }
    }
    //
    // Increment in fine grid velocity to velocity registers.
    //
    if (level > 0)
    {
        const Real mult = 1.0/( (double) parent->MaxRefRatio(level-1) );
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            mac_reg[level].FineAdd(u_mac[dir],area[level][dir],dir,0,0,1,mult);
        }

        Real sumreg = mac_reg[level].SumReg(0);

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            cout << "LEVEL "
                 << level
                 << " MACREG: FineAdd sum = "
                 << sumreg << endl;
        }
    }
    //
    // If finest level possible no need to keep phi for boundary conditions.
    //
    if (level == max_level)
    {
        delete mac_phi;
        mac_phi = 0;
    }
}

//
// Compute the corrective pressure used in the mac_sync.
//

void
MacProj::mac_sync_solve (int       level,
                         MultiFab* u_mac, 
                         Real      dt,
                         MultiFab* rho_half,
                         IntVect&  fine_ratio)
{
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "... mac_sync_solve at level " << level << NL;
    }
    assert(level < finest_level);

    const BoxArray& grids = LevelData[level].boxArray();
    const Geometry& geom  = parent->Geom(level);
    IntVect crse_ratio    = (level > 0) ? 
        parent->refRatio(level-1) : IntVect::TheZeroVector();
    const Real* dx             = geom.CellSize();
    const BoxArray& fine_boxes = LevelData[level+1].boxArray();
    //
    // Reusing storage here, since there should be no more need for the
    // values in mac_phi at this level and mac_sync_phi only need to last
    // into the call to mac_sync_compute.  Hope this works...  (LHH).
    //
    MultiFab* mac_sync_phi = &mac_phi_crse[level];
    //
    // Alloc and define RHS by doing a reflux-like operation in coarse
    // grid cells adjacent to fine grids.  The values in these
    // cells should be SUM{MR/VOL} where the sum is taken over
    // all edges of a cell that adjoin fine grids, MR = value in
    // MAC register, VOL = cell volume.  All other cells have a
    // value of zero (including crse cells under fine grids).
    //
    MultiFab Rhs(grids,1,0,Fab_allocate);
    Rhs.setVal(0.0);
    //
    // Reflux subtracts values at hi edge of coarse cell and
    // adds values at lo edge.  We want the opposite here so
    // set scale to -1 & alloc space for Rhs.
    //
    FluxRegister& mr = mac_reg[level+1];
    const Real scale = -1.0;
    mr.Reflux(Rhs,volume[level],scale,0,0,1,geom);

    for (int kf = 0, nfine = fine_boxes.length(); kf < nfine; kf++)
    {
        Box bf = ::coarsen(fine_boxes[kf],fine_ratio);

        for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
        {
            assert(grids[Rhsmfi.index()] == Rhsmfi.validbox());

            if (Rhsmfi.validbox().intersects(bf))
            {
                Rhsmfi().setVal(0.0,(Rhsmfi.validbox() & bf),0);
            }
        }
    }
    //
    // Remove constant null space component from the rhs of the solve
    // when appropriate (i.e. when the grids span the whole domain AND
    // the boundary conditions are Neumann for phi on all sides of the domain.)
    // Note that Rhs does not yet have the radial scaling in it for r-z
    // problems so we must do explicit volume-weighting here.
    //
    if (fix_mac_sync_rhs)
    {
        int all_neumann = 1;
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (phys_bc->lo()[dir] == Outflow || phys_bc->hi()[dir] == Outflow)
                all_neumann = 0;
        }

        long size = 0;
        for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
        {
            size += Rhsmfi().box().numPts();
        }
        ParallelDescriptor::ReduceLongSum(size);

        if (size == geom.Domain().numPts() && all_neumann == 1)
        {
            Real sum = 0.0;
            Real vol = 0.0;
            FArrayBox vol_wgted_rhs;
            for (MultiFabIterator Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
            {
                vol_wgted_rhs.resize(Rhsmfi().box());
                vol_wgted_rhs.copy(Rhsmfi());
                DependentMultiFabIterator Volmfi(Rhsmfi, volume[level]);
                vol_wgted_rhs.mult(Volmfi());
                sum += vol_wgted_rhs.sum(0,1);
                vol += Volmfi().sum(0,1);
            }
            ParallelDescriptor::ReduceRealSum(sum);
            ParallelDescriptor::ReduceRealSum(vol);

            const Real fix = sum / vol;
            if (ParallelDescriptor::IOProcessor())
            {
                cout << "Average correction on mac sync RHS = " << fix << NL;
            }
            Rhs.plus(-fix, 0);
        }
    }

    mac_sync_phi->setVal(0.0);
    //
    // store the Dirichlet boundary condition for mac_sync_phi in mac_bndry
    //
    MacBndry mac_bndry(grids,1,geom);
    const int src_comp = 0;
    const int dest_comp = 0;
    const int num_comp = 1;
    if (level == 0)
    {
        mac_bndry.setBndryValues(*mac_sync_phi,src_comp,dest_comp,num_comp,
                                 *phys_bc);
    }
    else
    {
        BoxArray crse_boxes(grids);
        crse_boxes.coarsen(crse_ratio);
        const int n_ghost = 0;
        BndryRegister crse_br(crse_boxes,n_ghost,1,1,num_comp);
        crse_br.setVal(0);
        mac_bndry.setBndryValues(crse_br,src_comp,*mac_sync_phi,src_comp,
                                 dest_comp,num_comp,crse_ratio, *phys_bc);
    }
    //
    // Now define edge centered coefficients and adjust RHS for MAC solve.
    //
    const Real rhs_scale = 2.0/dt;
    //
    // Solve the sync system.
    //
    mac_sync_driver(mac_bndry, grids, use_cg_solve, level, dx, dt,
                    mac_sync_tol, mac_abs_tol, rhs_scale, area[level],
                    volume[level], Rhs, rho_half, u_mac, mac_sync_phi);
}

//
// After solving for mac_sync_phi in mac_sync_solve(), we
// can now do the sync advect step.  This consists of two steps
//
// 1. compute u_corr as the gradient of mac_sync_phi
// 2. compute advective tendency of u_corr and
//    add into Vsync or Ssync
//
// If increment_sync is non-null, the (i-BL_SPACEDIM)-th component 
// of (*Ssync) is incremented only when increment[i]==1
// This is useful if that component gets incrmnted in a non-standard way.
//

void
MacProj::mac_sync_compute (int           level,
                           MultiFab*     u_mac, 
                           MultiFab*     Vsync,
                           MultiFab*     Ssync,
                           MultiFab*     rho_half,
                           FluxRegister* adv_flux_reg,
                           Array<int>    is_conservative,
                           Real          prev_time, 
                           Real          pres_prev_time,
                           Real          dt, 
                           int           NUM_STATE,
                           Real          be_cn_theta,
                           const int*    increment_sync)
{
    FArrayBox Rho, tforces, tvelforces;
    FArrayBox xflux, yflux, zflux, Gp;
    FArrayBox grad_phi[BL_SPACEDIM];
    //
    // Get parameters.
    //
    const BoxArray& grids    = LevelData[level].boxArray();
    const Geometry& geom     = parent->Geom(level);
    const Real* dx           = geom.CellSize();
    const int numscal        = NUM_STATE - BL_SPACEDIM;
    MultiFab* mac_sync_phi   = &mac_phi_crse[level];
    NavierStokes& ns_level   = *(NavierStokes*) &(parent->getLevel(level));
    Godunov* godunov         = ns_level.godunov;
    bool use_forces_in_trans = godunov->useForcesInTrans();

    MultiFab vel_visc_terms;

    if (use_forces_in_trans)
    {
        vel_visc_terms.define(grids,BL_SPACEDIM,1,Fab_allocate);

        if (be_cn_theta != 1.0)
            ns_level.getViscTerms(vel_visc_terms,Xvel,BL_SPACEDIM,prev_time);
        else
            vel_visc_terms.setVal(0,1);
    }
    //
    // Get viscous forcing.
    //
    MultiFab visc_terms(grids,NUM_STATE,1);

    if (be_cn_theta != 1.0)
        ns_level.getViscTerms(visc_terms,Xvel,NUM_STATE,prev_time);
    else
        visc_terms.setVal(0,1);

    Array<int> ns_level_bc, bndry[BL_SPACEDIM];
    //
    // FillPatch()d stuff allocated on heap ...
    //
    MultiFab* divu_fp = ns_level.getDivCond(1,prev_time);

    FillPatchIterator P_fpi(ns_level,ns_level.get_old_data(Press_Type),1,
                            pres_prev_time,Press_Type,0,1);

    FillPatchIterator S_fpi(ns_level,visc_terms,HYP_GROW,
                            prev_time,State_Type,0,NUM_STATE);
    //
    // Compute the mac sync correction.
    //
    for ( ; S_fpi.isValid() && P_fpi.isValid(); ++S_fpi, ++P_fpi)
    {
        DependentMultiFabIterator u_mac0mfi(S_fpi,u_mac[0]);
        DependentMultiFabIterator u_mac1mfi(S_fpi, u_mac[1]);
        DependentMultiFabIterator volumemfi(S_fpi,volume[level]);
        DependentMultiFabIterator area0mfi(S_fpi,area[level][0]);
        DependentMultiFabIterator area1mfi(S_fpi,area[level][1]);
#if (BL_SPACEDIM == 3)
        DependentMultiFabIterator area2mfi(S_fpi,area[level][2]);
        DependentMultiFabIterator u_mac2mfi(S_fpi,u_mac[2]);
#endif
        DependentMultiFabIterator Vsyncmfi(S_fpi,*Vsync);
        DependentMultiFabIterator Ssyncmfi(S_fpi,*Ssync);
        DependentMultiFabIterator rho_halfmfi(S_fpi,*rho_half);
        DependentMultiFabIterator mac_sync_phimfi(S_fpi,*mac_sync_phi);
        DependentMultiFabIterator visc_termsmfi(S_fpi,visc_terms);

        const int i     = S_fpi.index();
        FArrayBox& S    = S_fpi();
        FArrayBox& divu = (*divu_fp)[i];
        //
        // Step 1: compute ucorr = grad(phi)/rhonph
        //
        // Create storage for corrective velocities.
        //
        grad_phi[0].resize(::surroundingNodes(grids[i],0),1);
        grad_phi[1].resize(::surroundingNodes(grids[i],1),1);
#if (BL_SPACEDIM == 3)
        grad_phi[2].resize(::surroundingNodes(grids[i],2),1);
#endif
        mac_vel_update(1, grad_phi[0], grad_phi[1],
#if (BL_SPACEDIM == 3)
                       grad_phi[2],
#endif
                       mac_sync_phimfi(), &rho_halfmfi(),
                       0, grids[i], level, i, dx, dt/2.0);
        //
        // Step 2: compute Mac correction by calling GODUNOV box
        //
        // Get needed data.
        //
        Rho.resize(::grow(grids[i],1),1);
        Rho.copy(S,Density,0,1);

        ns_level.getForce(tforces,i,1,0,NUM_STATE,Rho);
        //
        // Compute total forcing terms.
        //
        ns_level.getGradP(P_fpi(),Gp,grids[i],1);

        godunov->Sum_tf_gp_visc(tforces, visc_termsmfi(), Gp, Rho);
        godunov->Sum_tf_divu_visc(S, tforces, BL_SPACEDIM, numscal,
                                  visc_termsmfi(), BL_SPACEDIM, divu, Rho, 1);

        if (use_forces_in_trans)
        {
            DependentMultiFabIterator dmfi(S_fpi, vel_visc_terms);
            ns_level.getForce(tvelforces,i,1,Xvel,BL_SPACEDIM,Rho);
            godunov->Sum_tf_gp_visc(tvelforces, dmfi(), Gp, Rho);
        }
        //
        // Set up the workspace for the godunov Box.
        //
        D_TERM(bndry[0] = ns_level.getBCArray(State_Type,i,0,1);,
               bndry[1] = ns_level.getBCArray(State_Type,i,1,1);,
               bndry[2] = ns_level.getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 0,
                       xflux, bndry[0].dataPtr(),
                       yflux, bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)
                       zflux, bndry[2].dataPtr(),
#endif
                       S, Rho, tvelforces);
        //
        // Get the sync FABS.
        //
        FArrayBox& u_sync = Vsyncmfi();
        FArrayBox& s_sync = Ssyncmfi();
        //
        // Loop over state components and compute the sync advective component.
        //
        for (int comp = 0; comp < NUM_STATE; comp++)
        {
            int do_comp = (increment_sync == NULL);
            if (!do_comp)
            {
                do_comp = comp < BL_SPACEDIM || (increment_sync[comp]==1);
            }
            if (!do_comp)
                continue;
            const int u_ind    = comp;
            const int s_ind    = comp-BL_SPACEDIM;
            const int sync_ind = (comp < BL_SPACEDIM ? u_ind  : s_ind);
            FArrayBox& temp    = (comp < BL_SPACEDIM ? u_sync : s_sync);
            ns_level_bc        = ns_level.getBCArray(State_Type,i,comp,1);

            godunov->SyncAdvect(grids[i], dx, dt, level, area0mfi(), u_mac0mfi(),
                                grad_phi[0], xflux, area1mfi(), u_mac1mfi(),
                                grad_phi[1], yflux,
#if (BL_SPACEDIM == 3)                            
                                area2mfi(), u_mac2mfi(), grad_phi[2], zflux,
#endif
                                S, tforces, comp, temp, sync_ind,
                                is_conservative[comp], comp,
                                ns_level_bc.dataPtr(), volumemfi());
            //
            // NOTE: the signs here are opposite from VELGOD.
            // NOTE: fluxes expected to be in extensive form.
            //
            if (level > 0)
            {
                adv_flux_reg->FineAdd(xflux,0,i,0,comp,1,-dt);
                adv_flux_reg->FineAdd(yflux,1,i,0,comp,1,-dt);
#if (BL_SPACEDIM == 3)
                adv_flux_reg->FineAdd(zflux,2,i,0,comp,1,-dt);
#endif
            }
        }
        //
        // Include grad_phi in the mac registers corresponding
        // to the next coarsest interface.
        //
        if (level > 0)
        {
            const Real mult =  -1.0/( (double) parent->MaxRefRatio(level-1));
            mac_reg[level].FineAdd(grad_phi[0],area0mfi(),0,i,0,0,1,mult);
            mac_reg[level].FineAdd(grad_phi[1],area1mfi(),1,i,0,0,1,mult);
#if (BL_SPACEDIM == 3)
            mac_reg[level].FineAdd(grad_phi[2],area2mfi(),2,i,0,0,1,mult);
#endif
        }
        //
        // Multiply the sync term by dt -- now done in the calling routine.
        //
    }

    delete divu_fp;
}

//
// This routine does a sync advect step for a single 
// scalar component. Unlike the preceding routine, the
// half-time edge states are passed in from the calling routine.
// This routine is useful when the edge states are computed
// in a physics-class-specific manner. (For example, as they are
// in the calculation of div rho U h = div U sum_l (rho Y)_l h_l(T)).
//

void
MacProj::mac_sync_compute (int           level,
                           MultiFab*     u_mac, 
                           MultiFab*     Ssync,
                           int           comp,
                           const MultiFab* const* sync_edges,
			   int           eComp,
                           MultiFab*     rho_half,
                           FluxRegister* adv_flux_reg,
                           Array<int>    is_conservative, 
                           Real          dt)
{
    assert(comp >= BL_SPACEDIM);

    FArrayBox xflux, yflux, zflux;
    FArrayBox grad_phi[BL_SPACEDIM];

    const int s_ind        = comp - BL_SPACEDIM;    
    const BoxArray& grids  = LevelData[level].boxArray();
    const Geometry& geom   = parent->Geom(level);
    MultiFab* mac_sync_phi = &mac_phi_crse[level];
    NavierStokes& ns_level = *(NavierStokes*) &(parent->getLevel(level));

    Godunov godunov(512);
    //
    // Compute the mac sync correction.
    //
    for (MultiFabIterator Ssyncmfi(*Ssync); Ssyncmfi.isValid(); ++Ssyncmfi)
    {
        DependentMultiFabIterator volumemfi(Ssyncmfi, volume[level]);
        DependentMultiFabIterator area0mfi(Ssyncmfi, area[level][0]);
        DependentMultiFabIterator area1mfi(Ssyncmfi, area[level][1]);
        DependentMultiFabIterator sync_edges0mfi(Ssyncmfi, *sync_edges[0]);
        DependentMultiFabIterator sync_edges1mfi(Ssyncmfi, *sync_edges[1]);
#if (BL_SPACEDIM == 3)
        DependentMultiFabIterator area2mfi(Ssyncmfi, area[level][2]);
        DependentMultiFabIterator sync_edges2mfi(Ssyncmfi, *sync_edges[2]);
#endif
        DependentMultiFabIterator rho_halfmfi(Ssyncmfi, *rho_half);
        DependentMultiFabIterator mac_sync_phimfi(Ssyncmfi, *mac_sync_phi);
        //
        // Step 1: compute ucorr = grad(phi)/rhonph
        //
        grad_phi[0].resize(::surroundingNodes(grids[Ssyncmfi.index()],0),1);
        grad_phi[1].resize(::surroundingNodes(grids[Ssyncmfi.index()],1),1);
#if (BL_SPACEDIM == 3)
        grad_phi[2].resize(::surroundingNodes(grids[Ssyncmfi.index()],2),1);
#endif
        mac_vel_update(1,
                       grad_phi[0],
                       grad_phi[1],
#if (BL_SPACEDIM == 3)
                       grad_phi[2],
#endif
                       mac_sync_phimfi(),
                       &rho_halfmfi(), 0,
                       grids[Ssyncmfi.index()], level, Ssyncmfi.index(),
                       geom.CellSize(), dt/2.0);
        //
        // Step 2: compute Mac correction by advecting the edge states.
        //
        xflux.resize(::surroundingNodes(grids[Ssyncmfi.index()],0),1);
        yflux.resize(::surroundingNodes(grids[Ssyncmfi.index()],1),1);
        xflux.copy(sync_edges0mfi(),eComp,0,1);
        yflux.copy(sync_edges1mfi(),eComp,0,1);
#if (BL_SPACEDIM == 3)
        zflux.resize(::surroundingNodes(grids[Ssyncmfi.index()],2),1);
        zflux.copy(sync_edges2mfi(),eComp,0,1);
#endif
        godunov.ComputeSyncAofs(grids[Ssyncmfi.index()],
                                area0mfi(),
                                grad_phi[0],       xflux,
                                
                                area1mfi(),
                                grad_phi[1],       yflux,
#if (BL_SPACEDIM == 3)                            
                                area2mfi(),
                                grad_phi[2],       zflux,
#endif
                                volumemfi(), Ssyncmfi(), s_ind, 
                                is_conservative[comp]);
        //
        // NOTE: the signs here are opposite from VELGOD.
        // NOTE: fluxes expected to be in extensive form.
        //
        if (level > 0)
        {
            adv_flux_reg->FineAdd(xflux,0,Ssyncmfi.index(),0,comp,1,-dt);
            adv_flux_reg->FineAdd(yflux,1,Ssyncmfi.index(),0,comp,1,-dt);
#if (BL_SPACEDIM == 3)
            adv_flux_reg->FineAdd(zflux,2,Ssyncmfi.index(),0,comp,1,-dt);
#endif
        }
        //
        // Multiply the sync term by dt -- now done in the calling routine.
        //
    }
}

//
// Check the mac divergence.
//

void
MacProj::check_div_cond (int      level,
                         MultiFab U_edge[]) const
{
    const BoxArray& grids = LevelData[level].boxArray();

    Real sum = 0.0;

    FArrayBox dmac;

    for (MultiFabIterator U_edge0mfi(U_edge[0]); U_edge0mfi.isValid();
         ++U_edge0mfi)
    {
        DependentMultiFabIterator U_edge1mfi(U_edge0mfi, U_edge[1]);
        DependentMultiFabIterator volumemfi(U_edge0mfi, volume[level]);
        DependentMultiFabIterator area0mfi(U_edge0mfi, area[level][0]);
        DependentMultiFabIterator area1mfi(U_edge0mfi, area[level][1]);

#if (BL_SPACEDIM == 3)
        DependentMultiFabIterator U_edge2mfi(U_edge0mfi, U_edge[2]);
        DependentMultiFabIterator area2mfi(U_edge0mfi, area[level][2]);
#endif
        dmac.resize(grids[U_edge0mfi.index()],1);

        const FArrayBox& uxedge = U_edge0mfi();
        const FArrayBox& uyedge = U_edge1mfi();
        const FArrayBox& xarea  = area0mfi();
        const FArrayBox& yarea  = area1mfi();
        const FArrayBox& vol    = volumemfi();

        DEF_LIMITS(dmac,dmac_dat,dlo,dhi);
        DEF_CLIMITS(uxedge,ux_dat,uxlo,uxhi);
        DEF_CLIMITS(uyedge,uy_dat,uylo,uyhi);
        DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
        DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);
        DEF_CLIMITS(vol,vol_dat,vlo,vhi);

#if (BL_SPACEDIM == 2)
        FORT_MACDIV(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
                    ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi), 
                    ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                    vol_dat,ARLIM(vlo),ARLIM(vhi));
#endif

#if (BL_SPACEDIM == 3)
        const FArrayBox& uzedge = U_edge2mfi();
        DEF_CLIMITS(uzedge,uz_dat,uzlo,uzhi);
        const FArrayBox& zarea = area2mfi();
        DEF_CLIMITS(zarea,az_dat,azlo,azhi);

        FORT_MACDIV(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
                    ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
                    uz_dat,ARLIM(uzlo),ARLIM(uzhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi),
                    ay_dat,ARLIM(aylo),ARLIM(ayhi),
                    az_dat,ARLIM(azlo),ARLIM(azhi),
                    vol_dat,ARLIM(vlo),ARLIM(vhi));
#endif
        sum += dmac.sum(0);
    }

    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    ParallelDescriptor::ReduceRealSum(sum,IOProc);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "SUM of DIV(U_edge) = " << sum << NL;
    }
}

void
MacProj::set_outflow_bcs (int             level,
                          MultiFab*       mac_phi,
                          const MultiFab* u_mac, 
                          const MultiFab& S,
                          const MultiFab& divu)
{
    //
    // Warning: This code looks about right, but hasn't really been tested yet.
    //
    // This code is very similar to the outflow BC stuff in the Projection
    // class except that here the the phi to be solved for lives on the
    // out-directed faces.  The projection equation to satisfy is
    //
    //   (1/r)(d/dr)[r/rho dphi/dr] = dv/dr - S
    //
    // Outflow here is assumed to occur only at yhi faces (should be easy to
    // generalize, however).
    //
#if (BL_SPACEDIM == 2)
    const Orientation outFace(1, Orientation::high);

    const int rzflag  = CoordSys::IsRZ();
    const Real* dx    = parent->Geom(level).CellSize();
    const Box& domain = parent->Geom(level).Domain();
    //
    // Load cc data (rho, divu).
    //
    const int bndBxWdth = 1;
    const int dir       = outFace.coordDir();

    Box ccBndBox = ::adjCell(domain,outFace,bndBxWdth).shift(dir,-bndBxWdth);
    Box phiBox   = ::adjCell(domain,outFace,1);
    //
    // Note: this code assumes that we can fill ccBndBox from S.
    //
    // TODO -- how to generalize this ???
    //
    assert(::complementIn(ccBndBox,BoxList(S.boxArray())).isEmpty());

    FArrayBox rhodat(ccBndBox,1);
    FArrayBox divudat(ccBndBox,1);
    FArrayBox phidat(phiBox,1);

    const BoxArray& ba = S.boxArray();

    assert(ba == divu.boxArray());
    //
    // Fill rhodat & divudat.
    //
    S.copy(rhodat, Density, 0, 1);
    divu.copy(divudat, 0, 0, 1);
    //
    // Load ec data.
    //
    FArrayBox uedat[BL_SPACEDIM-1];

    for (int i = 0, cnt = 0; i < BL_SPACEDIM; ++i)
	if (i != dir)
	    uedat[cnt++].resize(::surroundingNodes(ccBndBox,i), 1);

    for (int i = 0, cnt = 0; i < BL_SPACEDIM; ++i)
	if (i != dir)
            u_mac[i].copy(uedat[cnt++]);
    //
    // Make cc r (set = 1 if cartesian).
    //
    Array<Real> rcen(ccBndBox.length(0), 1.0);

    if (CoordSys::IsRZ() == 1) 
	parent->Geom(level).GetCellLoc(rcen, ccBndBox, 0);
    //
    // Compute boundary solution.
    //
    const int isPeriodicInX = parent->Geom(level).isPeriodic(0);

    FORT_MACPHIBC(ARLIM(uedat[0].loVect()),ARLIM(uedat[0].hiVect()),uedat[0].dataPtr(),
		  ARLIM(divudat.loVect()), ARLIM(divudat.hiVect()),divudat.dataPtr(),
		  ARLIM(rhodat.loVect()),  ARLIM(rhodat.hiVect()), rhodat.dataPtr(),
		  ARLIM(ccBndBox.loVect()),ARLIM(ccBndBox.hiVect()),
		  rcen.dataPtr(), &dx[0],
		  ARLIM(phidat.loVect()),  ARLIM(phidat.hiVect()), phidat.dataPtr(),
		  &isPeriodicInX);

    for (MultiFabIterator mfi(*mac_phi); mfi.isValid(); ++mfi)
    {
        if (mfi.validbox().intersects(phidat.box()))
        {
            Box ovlp = mfi.validbox() & phidat.box();
	    mfi().copy(phidat,ovlp);
        }
    }
#else
    BoxLib::Error("MacProj::set_outflow_bcs(): not implemented yet for 3D");
#endif
}
