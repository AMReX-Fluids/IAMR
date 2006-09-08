
//
// $Id: MacProj.cpp,v 1.105 2006-09-08 21:21:02 almgren Exp $
//
#include <winstd.H>

#include <LO_BCTYPES.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <MacProj.H>
#include <MacBndry.H>
#include <MacOpMacDrivers.H>
#include <PorousMedia.H>
#include <MACPROJ_F.H>
#include <MacOutFlowBC.H>

#ifndef _PorousMedia_H_
enum StateType {State_Type=0, Press_Type};
enum StateNames  { N_w, N_o, Enthalpy};
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define DEF_BOX_LIMITS(box,boxlo,boxhi)      \
const int* boxlo = (box).loVect();           \
const int* boxhi = (box).hiVect();

#define GEOM_GROW 1
#define HYP_GROW 3

int  MacProj::verbose          = 0;
bool MacProj::use_cg_solve     = false;
namespace
{
  bool        use_hypre_solve  = false;
  bool        use_fboxlib_mg   = false;
}
Real MacProj::mac_tol          = 1.0e-12;
Real MacProj::mac_abs_tol      = 1.0e-16;
Real MacProj::mac_sync_tol     = 1.0e-8;
int  MacProj::do_outflow_bcs   = 1;
int  MacProj::fix_mac_sync_rhs = 0;
int  MacProj::check_umac_periodicity = 1;

namespace
{
  Real umac_periodic_test_Tol    = 1.e-10;
}

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
    anel_coeff(_finest_level+1), radius(_finest_level+1) 
{
    read_params();

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Creating mac_projector\n";

    finest_level_allocated = finest_level;

    for (int lev = 0; lev <= finest_level; lev++)
       anel_coeff[lev] = 0;
}

MacProj::~MacProj () {}

void
MacProj::read_params ()
{
    ParmParse pp("mac");

    pp.query( "v",                verbose          );
    pp.query( "mac_tol",          mac_tol          );
    pp.query( "mac_sync_tol",     mac_sync_tol     );
    pp.query( "use_cg_solve",     use_cg_solve     );
#if MG_USE_HYPRE
    pp.query( "use_hypre_solve",     use_hypre_solve);
#endif
#if MG_USE_FBOXLIB
    pp.query( "use_fboxlib_mg",     use_fboxlib_mg);
#endif
    pp.query( "mac_abs_tol",      mac_abs_tol      );
    pp.query( "do_outflow_bcs",   do_outflow_bcs   );
    pp.query( "fix_mac_sync_rhs", fix_mac_sync_rhs );
    pp.query("check_umac_periodicity",check_umac_periodicity);
    pp.query("umac_periodic_test_Tol",   umac_periodic_test_Tol);

    if ( use_cg_solve && use_hypre_solve )
      {
	BoxLib::Error("MacProj::read_params: cg_solve && .not. hypre_solve");
      }
    if ( use_cg_solve && use_fboxlib_mg )
      {
	BoxLib::Error("MacProj::read_params: cg_solve && .not. fboxlib_solve");
      }
}

void
MacProj::install_level (int                   level,
                        AmrLevel*             level_data,
                        MultiFab&             _volume,
                        MultiFab*             _area,
                        Array< Array<Real> >* _radius)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Installing MacProj level " << level << '\n';

    if (parent->finestLevel() < finest_level)
        for (int lev = parent->finestLevel() + 1; lev <= finest_level; lev++)
            mac_reg.clear(lev);

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
        mac_reg.clear(level);
        mac_reg.set(level,new FluxRegister(LevelData[level].boxArray(),
                                           parent->refRatio(level-1),level,1));
    }

    if (level > anel_coeff.size()-1) {
       anel_coeff.resize(level+1);
       anel_coeff[level] = 0;
    }
}
void
MacProj::install_anelastic_coefficient (int               level,
                                        Real**            _anel_coeff)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Installing anel_coeff into MacProj level " << level << '\n';

    if (level > anel_coeff.size()-1) anel_coeff.resize(level+1);
    anel_coeff.set(level, _anel_coeff);
}

void
MacProj::BuildPhiBC (int level)
{
    const BoxArray& grids   = LevelData[level].boxArray();
    const Geometry& geom    = parent->Geom(level);
    const int       ngrds   = grids.size();
    phi_bcs[level].resize(ngrds);
    const Box&      domain  = geom.Domain();
    const int*      domlo   = domain.loVect();
    const int*      domhi   = domain.hiVect();
    const int*      phys_lo = phys_bc->lo();
    const int*      phys_hi = phys_bc->hi();

    for (int i = 0; i < ngrds; i++)
    {
        BCRec&     bc = phi_bcs[level][i];
        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (lo[dir] == domlo[dir])
            {
                bc.setLo(dir,phys_lo[dir]==Outflow ? LO_DIRICHLET : LO_NEUMANN);
            }
            else
            {
                bc.setLo(dir,LO_DIRICHLET);
            }
            if (hi[dir] == domhi[dir])
            {
                bc.setHi(dir,phys_hi[dir]==Outflow ? LO_DIRICHLET : LO_NEUMANN);
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
            mac_phi_crse.set(level,new MultiFab(grids,1,1));
            mac_phi_crse[level].setVal(0.0);
        }
    }
}

void
MacProj::cleanup (int level)
{
    if (level < parent->maxLevel())
        mac_phi_crse.clear(level);
}

//
// Projection functions follow ...
//
static
bool
grids_on_side_of_domain (const BoxArray&    grids,
                         const Box&         domain,
                         const Orientation& outFace)
{
    const int idir = outFace.coordDir();

    if (outFace.isLow())
    {
        for (int igrid = 0; igrid < grids.size(); igrid++)
        { 
            if (grids[igrid].smallEnd(idir) == domain.smallEnd(idir))
            { 
                return true;
            }
        }
    }
  
    if (outFace.isHigh())
    {
        for (int igrid = 0; igrid < grids.size(); igrid++)
        {
            if (grids[igrid].bigEnd(idir) == domain.bigEnd(idir))
            {
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
                      MultiFab&       mac_coef,
                      Real            dt,
                      Real            time,
                      const MultiFab& divu,
                      int             have_divu)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::mac_project()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mac_project at level " << level << '\n';

    const BoxArray& grids      = LevelData[level].boxArray();
    const Geometry& geom       = parent->Geom(level);
    const Real*     dx         = geom.CellSize();
    const int       max_level  = parent->maxLevel();
    MultiFab*       mac_phi    = 0;
    PorousMedia&    ns        = *(PorousMedia*) &(parent->getLevel(level));
    IntVect         crse_ratio = level > 0 ? parent->refRatio(level-1)
                                           : IntVect::TheZeroVector();
    //
    // If finest level possible no need to make permanent mac_phi for bcs.
    //
    if (level == max_level)
        mac_phi = new MultiFab(grids,1,1);
    else
        mac_phi = &mac_phi_crse[level];

    mac_phi->setVal(0.0);
    //
    // HACK!!!
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
        //const int extent_rad = 1;
        const int extent_rad = 2;
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
    MultiFab Rhs(grids,1,0);
    Rhs.copy(divu);

    int the_solver = 0;
    if (use_cg_solve)
    {
	the_solver = 1;
    }
    else if ( use_fboxlib_mg )
    {
	the_solver = 3;
    }

    if (anel_coeff[level] != 0) scaleArea(level,area[level],anel_coeff[level]);

    mac_level_driver(mac_bndry, *phys_bc, grids, the_solver, level, 
                     dx, dt, mac_tol, mac_abs_tol, rhs_scale, 
                     area[level], volume[level], mac_coef, Rhs, u_mac, mac_phi);
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

        if (verbose)
        {
            Real sumreg =  mr.SumReg(0);

            if (ParallelDescriptor::IOProcessor())
            {
                std::cout << "LEVEL "                   << level
                          << " MACREG: CrseInit sum = " << sumreg << std::endl;
            }
        }
    }
    //
    // Increment in fine grid velocity to velocity registers.
    //
    if (level > 0)
    {
        const Real mult = 1.0/parent->nCycle(level);

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            mac_reg[level].FineAdd(u_mac[dir],area[level][dir],dir,0,0,1,mult);
        }

        if (verbose)
        {
            Real sumreg = mac_reg[level].SumReg(0);

            if (ParallelDescriptor::IOProcessor())
            {
                std::cout << "LEVEL "                  << level
                          << " MACREG: FineAdd sum = " << sumreg << std::endl;
            }
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

    if (check_umac_periodicity)
        test_umac_periodic(level,u_mac);

    if (anel_coeff[level] != 0) rescaleArea(level,area[level],anel_coeff[level]);
}

//
// Compute the corrective pressure used in the mac_sync.
//

void
MacProj::mac_sync_solve (int       level,
                         Real      dt,
                         MultiFab* mac_coef,
                         IntVect&  fine_ratio)
{
    BL_ASSERT(level < finest_level);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mac_sync_solve at level " << level << '\n';

    const BoxArray& grids      = LevelData[level].boxArray();
    const Geometry& geom       = parent->Geom(level);
    const Real*     dx         = geom.CellSize();
    const BoxArray& fine_boxes = LevelData[level+1].boxArray();
    IntVect         crse_ratio = level > 0 ? parent->refRatio(level-1)
                                           : IntVect::TheZeroVector();
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
    MultiFab Rhs(grids,1,0);
    Rhs.setVal(0.0);
    //
    // Reflux subtracts values at hi edge of coarse cell and
    // adds values at lo edge.  We want the opposite here so
    // set scale to -1 & alloc space for Rhs.
    //
    FluxRegister& mr = mac_reg[level+1];
    const Real scale = -1.0;
    mr.Reflux(Rhs,volume[level],scale,0,0,1,geom);

    for (int kf = 0, nfine = fine_boxes.size(); kf < nfine; kf++)
    {
        Box bf = BoxLib::coarsen(fine_boxes[kf],fine_ratio);

        for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
        {
            BL_ASSERT(grids[Rhsmfi.index()] == Rhsmfi.validbox());

            Box isect = Rhsmfi.validbox() & bf;

            if (isect.ok())
            {
                Rhs[Rhsmfi].setVal(0.0,isect,0);
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

        if (Rhs.boxArray().contains(geom.Domain()) && all_neumann == 1)
        {
            Real sum = 0.0;
            Real vol = 0.0;
            FArrayBox vol_wgted_rhs;
            for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
            {
                vol_wgted_rhs.resize(Rhs[Rhsmfi].box());
                vol_wgted_rhs.copy(Rhs[Rhsmfi]);
                vol_wgted_rhs.mult(volume[level][Rhsmfi]);
                sum += vol_wgted_rhs.sum(0,1);
                vol += volume[level][Rhsmfi].sum(Rhs[Rhsmfi].box(),0,1);
            }
            ParallelDescriptor::ReduceRealSum(sum);
            ParallelDescriptor::ReduceRealSum(vol);

            const Real fix = sum / vol;

            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << "Average correction on mac sync RHS = " << fix << '\n';

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
        const int in_rad     = 0;
        const int out_rad    = 1;
        //const int extent_rad = 1;
        const int extent_rad = 2;
        BndryRegister crse_br(crse_boxes,in_rad,out_rad,extent_rad,num_comp);
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
    int the_solver = 0;
    if (use_cg_solve)
    {
	the_solver = 1;
    }
    else if ( use_fboxlib_mg )
    {
	the_solver = 3;
    }
    if (anel_coeff[level] != 0) scaleArea(level,area[level],anel_coeff[level]);
    mac_sync_driver(mac_bndry, *phys_bc, grids, the_solver, level, dx, dt,
                    mac_sync_tol, mac_abs_tol, rhs_scale, area[level],
                    volume[level], Rhs, mac_coef, mac_sync_phi);
    if (anel_coeff[level] != 0) rescaleArea(level,area[level],anel_coeff[level]);
}

//
// After solving for mac_sync_phi in mac_sync_solve(), we
// can now do the sync advect step.  This consists of two steps
//
// 1. compute u_corr as the gradient of mac_sync_phi
// 2. compute advective tendency of u_corr and
//    add into Ssync
//
// If increment_sync is non-null, the (i-BL_SPACEDIM)-th component 
// of (*Ssync) is incremented only when increment[i]==1
// This is useful if that component gets incrmnted in a non-standard way.
//

void
MacProj::mac_sync_compute (int                   level,
                           MultiFab*             u_mac, 
                           MultiFab*             Ssync,
                           MultiFab*             mac_coef,
                           FluxRegister*         adv_flux_reg,
                           Array<AdvectionForm>& advectionType,
                           Real                  prev_time, 
                           Real                  prev_pres_time,
                           Real                  dt, 
                           int                   NUM_STATE,
                           Real                  be_cn_theta,
                           const Array<int>&     increment_sync)
{
    FArrayBox tforces;
    FArrayBox xflux, yflux, zflux;
    FArrayBox grad_phi[BL_SPACEDIM];
    //
    // Get parameters.
    //
    const BoxArray& grids               = LevelData[level].boxArray();
    const Geometry& geom                = parent->Geom(level);
    const Real*     dx                  = geom.CellSize();
    const int       numscal             = NUM_STATE - BL_SPACEDIM;
    MultiFab*       mac_sync_phi        = &mac_phi_crse[level];
    PorousMedia&   ns_level            = *(PorousMedia*) &(parent->getLevel(level));
    Godunov*        godunov             = ns_level.godunov;
    bool            use_forces_in_trans = godunov->useForcesInTrans()?true:false;

    MultiFab scal_visc_terms(grids,numscal,1);

    scal_visc_terms.setVal(0,1); // Initialize to make calls below safe
    //
    // Get viscous forcing.
    //
    if (be_cn_theta != 1.0) 
    {
        int i;
        bool do_get_visc_terms = false;

        for (i=0; i < BL_SPACEDIM; ++i)
            if (!increment_sync.size() || increment_sync[i]==1)
                do_get_visc_terms = true;

        do_get_visc_terms = false;
        for (i=BL_SPACEDIM; i < increment_sync.size(); ++i)
            if (!increment_sync.size() || increment_sync[i]==1)
                do_get_visc_terms = true;

        if (do_get_visc_terms)
            ns_level.getViscTerms(scal_visc_terms,BL_SPACEDIM,numscal,prev_time);
    }

    Array<int> ns_level_bc, bndry[BL_SPACEDIM];
    //
    // FillPatch()d stuff allocated on heap ...
    //
    MultiFab Gp(grids,BL_SPACEDIM,1);

    ns_level.getGradP(Gp, prev_pres_time);

    MultiFab* divu_fp = ns_level.getDivCond(1,prev_time);

    FluxRegister* temp_reg = 0;

    //
    // Compute the mac sync correction.
    //
    for (FillPatchIterator P_fpi(ns_level,ns_level.get_old_data(Press_Type),1,prev_pres_time,Press_Type,0,1),
                           S_fpi(ns_level,scal_visc_terms,HYP_GROW,prev_time,State_Type,0,NUM_STATE);
         S_fpi.isValid() && P_fpi.isValid();
         ++S_fpi, ++P_fpi)
    {
        const int i     = S_fpi.index();
        FArrayBox& S    = S_fpi();
        FArrayBox& divu = (*divu_fp)[i];

        FArrayBox U;

        U.resize(S.box(),BL_SPACEDIM);
        U.copy(S_fpi(),0,0,BL_SPACEDIM);
        //
        // Step 1: compute ucorr = grad(phi)/rhonph
        //
        // Create storage for corrective velocities.
        //
        D_TERM(grad_phi[0].resize(BoxLib::surroundingNodes(grids[i],0),1);,
               grad_phi[1].resize(BoxLib::surroundingNodes(grids[i],1),1);,
               grad_phi[2].resize(BoxLib::surroundingNodes(grids[i],2),1););

        mac_vel_update(1,D_DECL(grad_phi[0],grad_phi[1],grad_phi[2]),
                       (*mac_sync_phi)[S_fpi], &(*mac_coef)[S_fpi],
                       grids[i], level, i, dx, dt/2.0);
        //
        // Step 2: compute Mac correction by calling GODUNOV box
        ns_level.getForce(tforces,i,1,0,NUM_STATE);
        //
        // Compute total forcing terms.
        //
        godunov->Sum_tf_divu_visc(S, tforces, 0, numscal,
                                  scal_visc_terms[S_fpi], 0, divu, 1);

        //
        // Set up the workspace for the godunov Box.
        //
        D_TERM(bndry[0] = ns_level.getBCArray(State_Type,i,0,1);,
               bndry[1] = ns_level.getBCArray(State_Type,i,1,1);,
               bndry[2] = ns_level.getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 0,
                       xflux, bndry[0].dataPtr(),
#if (BL_SPACEDIM == 2)
                       yflux, bndry[1].dataPtr());
#elif (BL_SPACEDIM == 3)
                       yflux, bndry[1].dataPtr(),
                       zflux, bndry[2].dataPtr());
#endif
        //
        // Get the sync FABS.
        //
        FArrayBox& s_sync = (*Ssync)[S_fpi];
        //
        // Loop over state components and compute the sync advective component.
        //
        for (int comp = 0; comp < NUM_STATE; comp++)
        {
            if (!increment_sync.size() || increment_sync[comp]==1)
            {
                const int  sync_ind = comp;
                FArrayBox& temp     = s_sync;
                ns_level_bc         = ns_level.getBCArray(State_Type,i,comp,1);

                int use_conserv_diff = (advectionType[comp] == Conservative) ? true : false;

                godunov->SyncAdvect(grids[i], dx, dt, level, 
                                    area[level][0][S_fpi], u_mac[0][S_fpi], grad_phi[0], xflux, 
                                    area[level][1][S_fpi], u_mac[1][S_fpi], grad_phi[1], yflux,
#if (BL_SPACEDIM == 3)                            
                                    area[level][2][S_fpi], u_mac[2][S_fpi], grad_phi[2], zflux,
#endif
                                    S, tforces, divu, comp, temp, sync_ind,
                                    use_conserv_diff, comp,
                                    ns_level_bc.dataPtr(), PRE_MAC, volume[level][S_fpi]);
                //
                // NOTE: the signs here are opposite from VELGOD.
                // NOTE: fluxes expected to be in extensive form.
                //
                if (level > 0)
                {
                    D_TERM(adv_flux_reg->FineAdd(xflux,0,i,0,comp,1,-dt);,
                           adv_flux_reg->FineAdd(yflux,1,i,0,comp,1,-dt);,
                           adv_flux_reg->FineAdd(zflux,2,i,0,comp,1,-dt););
                }
            }
        }
        //
        // Fill temp_reg with the normal fluxes.
        //
        int velpred = 0;
        //
        // Include grad_phi in the mac registers corresponding
        // to the next coarsest interface.
        //
        if (level > 0)
        {
            const Real mlt =  -1.0/( (double) parent->nCycle(level));
            D_TERM(mac_reg[level].FineAdd(grad_phi[0],area[level][0][S_fpi],0,i,0,0,1,mlt);,
                   mac_reg[level].FineAdd(grad_phi[1],area[level][1][S_fpi],1,i,0,0,1,mlt);,
                   mac_reg[level].FineAdd(grad_phi[2],area[level][2][S_fpi],2,i,0,0,1,mlt););
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
MacProj::mac_sync_compute (int                    level,
                           MultiFab*              Sync,
                           int                    comp,
                           int                    s_ind,
                           const MultiFab* const* sync_edges,
			   int                    eComp,
                           MultiFab*              mac_coef,
                           FluxRegister*          adv_flux_reg,
                           Array<AdvectionForm>&  advectionType, 
                           Real                   dt)
{
    FArrayBox xflux, yflux, zflux;
    FArrayBox grad_phi[BL_SPACEDIM];

    const BoxArray& grids        = LevelData[level].boxArray();
    const Geometry& geom         = parent->Geom(level);
    MultiFab*       mac_sync_phi = &mac_phi_crse[level];
    PorousMedia&   ns_level     = *(PorousMedia*) &(parent->getLevel(level));

    Godunov godunov(512);

    FluxRegister* temp_reg = 0;
    //
    // Compute the mac sync correction.
    //
    for (MFIter Syncmfi(*Sync); Syncmfi.isValid(); ++Syncmfi)
    {
        //
        // Step 1: compute ucorr = grad(phi)/rhonph
        //
        D_TERM(grad_phi[0].resize(BoxLib::surroundingNodes(grids[Syncmfi.index()],0),1);,
               grad_phi[1].resize(BoxLib::surroundingNodes(grids[Syncmfi.index()],1),1);,
               grad_phi[2].resize(BoxLib::surroundingNodes(grids[Syncmfi.index()],2),1););

        mac_vel_update(1,
                       D_DECL(grad_phi[0],grad_phi[1],grad_phi[2]),
                       (*mac_sync_phi)[Syncmfi],
                       &(*mac_coef)[Syncmfi], 
                       grids[Syncmfi.index()], level, Syncmfi.index(),
                       geom.CellSize(), dt/2.0);
        //
        // Step 2: compute Mac correction by advecting the edge states.
        //
        D_TERM(xflux.resize(BoxLib::surroundingNodes(grids[Syncmfi.index()],0),1);,
               yflux.resize(BoxLib::surroundingNodes(grids[Syncmfi.index()],1),1);,
               zflux.resize(BoxLib::surroundingNodes(grids[Syncmfi.index()],2),1););

        D_TERM(xflux.copy((*sync_edges[0])[Syncmfi],eComp,0,1);,
               yflux.copy((*sync_edges[1])[Syncmfi],eComp,0,1);,
               zflux.copy((*sync_edges[2])[Syncmfi],eComp,0,1););

        int use_conserv_diff = (advectionType[comp] == Conservative)
                                                             ? true : false;
        godunov.ComputeSyncAofs(grids[Syncmfi.index()],
                                area[level][0][Syncmfi],
                                grad_phi[0],       xflux,
                                
                                area[level][1][Syncmfi],
                                grad_phi[1],       yflux,
#if (BL_SPACEDIM == 3)                            
                                area[level][2][Syncmfi],
                                grad_phi[2],       zflux,
#endif
                                volume[level][Syncmfi], (*Sync)[Syncmfi],
                                s_ind, use_conserv_diff);
        //
        // NOTE: the signs here are opposite from VELGOD.
        // NOTE: fluxes expected to be in extensive form.
        //
        if (level > 0)
        {
            D_TERM(adv_flux_reg->FineAdd(xflux,0,Syncmfi.index(),0,comp,1,-dt);,
                   adv_flux_reg->FineAdd(yflux,1,Syncmfi.index(),0,comp,1,-dt);,
                   adv_flux_reg->FineAdd(zflux,2,Syncmfi.index(),0,comp,1,-dt););
        }
        //
        // Fill temp_reg with the normal fluxes.
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

    for (MFIter U_edge0mfi(U_edge[0]); U_edge0mfi.isValid(); ++U_edge0mfi)
    {
        dmac.resize(grids[U_edge0mfi.index()],1);

        const FArrayBox& uxedge = U_edge[0][U_edge0mfi];
        const FArrayBox& uyedge = U_edge[1][U_edge0mfi];
        const FArrayBox& xarea  = area[level][0][U_edge0mfi];
        const FArrayBox& yarea  = area[level][1][U_edge0mfi];
        const FArrayBox& vol    = volume[level][U_edge0mfi];

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
        const FArrayBox& uzedge = U_edge[2][U_edge0mfi];
        DEF_CLIMITS(uzedge,uz_dat,uzlo,uzhi);
        const FArrayBox& zarea = area[level][2][U_edge0mfi];
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

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealSum(sum,IOProc);
        
        if (ParallelDescriptor::IOProcessor())
            std::cout << "SUM of DIV(U_edge) = " << sum << '\n';
    }
}

//
// Structure used by test_umac_periodic().
//

struct TURec
{
    TURec ()
        :
        m_idx(-1),
        m_dim(-1)
    {}

    TURec (int        idx,
           int        dim,
           const Box& srcBox,
           const Box& dstBox)
        :
        m_srcBox(srcBox),
        m_dstBox(dstBox),
        m_idx(idx),
        m_dim(dim)
    {}

    FillBoxId m_fbid;
    Box       m_srcBox;
    Box       m_dstBox;
    int       m_idx;
    int       m_dim;
};

//
// Test that edge-based values agree across periodic boundary.
//

void
MacProj::test_umac_periodic (int level,MultiFab* u_mac)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::test_umac_periodic()");

    const Geometry& geom = parent->Geom(level);

    if (!geom.isAnyPeriodic()) return;

    const BoxArray&         grids  = LevelData[level].boxArray();
    const int               MyProc = ParallelDescriptor::MyProc();
    FArrayBox               diff;
    Array<IntVect>          pshifts(27);
    MultiFabCopyDescriptor  mfcd;
    std::vector<TURec>      pirm;
    MultiFabId              mfid[BL_SPACEDIM];

    for (int dim = 0; dim < BL_SPACEDIM; dim++)
    {
        if (geom.isPeriodic(dim))
        {
            Box eDomain = BoxLib::surroundingNodes(geom.Domain(),dim);

            mfid[dim] = mfcd.RegisterMultiFab(&u_mac[dim]);

            for (MFIter mfi(u_mac[dim]); mfi.isValid(); ++mfi)
            {
                Box eBox = u_mac[dim].boxArray()[mfi.index()];

                geom.periodicShift(eDomain, eBox, pshifts);

                for (int iiv = 0; iiv < pshifts.size(); iiv++)
                {
                    eBox += pshifts[iiv];

                    std::vector< std::pair<int,Box> > isects = u_mac[dim].boxArray().intersections(eBox);

                    for (int i = 0; i < isects.size(); i++)
                    {
                        const int j      = isects[i].first;
                        const Box srcBox = isects[i].second;
                        const Box dstBox = srcBox - pshifts[iiv];
                        TURec r(mfi.index(),dim,srcBox,dstBox);
                        r.m_fbid = mfcd.AddBox(mfid[dim],srcBox,0,j,0,0,1);
                        pirm.push_back(r);
                    }

                    eBox -= pshifts[iiv];
                }
            }
        }
    }

    mfcd.CollectData();

    for (int i = 0; i < pirm.size(); i++)
    {
        const int dim = pirm[i].m_dim;

        BL_ASSERT(pirm[i].m_fbid.box() == pirm[i].m_srcBox);
        BL_ASSERT(pirm[i].m_srcBox.sameSize(pirm[i].m_dstBox));
        BL_ASSERT(u_mac[dim].DistributionMap()[pirm[i].m_idx] == MyProc);

        diff.resize(pirm[i].m_srcBox, 1);

        mfcd.FillFab(mfid[dim], pirm[i].m_fbid, diff);

        diff.minus(u_mac[dim][pirm[i].m_idx],pirm[i].m_dstBox,diff.box(),0,0,1);

        const Real max_norm = diff.norm(0);

        if (max_norm > umac_periodic_test_Tol )
        {
            std::cout << "dir = "         << dim
                      << ", diff norm = " << max_norm
                      << " for region: "  << pirm[i].m_dstBox << std::endl;
            BoxLib::Error("Periodic bust in u_mac");
        }
    }
}

void
MacProj::scaleArea (int level, MultiFab* area, Real** anel_coeff)
{
    const BoxArray& grids = LevelData[level].boxArray();

    int mult = 1;

    for (MFIter mfi(*area); mfi.isValid(); ++mfi)
    {
        const FArrayBox& xarea  = area[0][mfi];
        const FArrayBox& yarea  = area[1][mfi];

        DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
        DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);

        const int* lo = grids[mfi.index()].loVect();
        const int* hi = grids[mfi.index()].hiVect();

        int anel_coeff_lo = lo[BL_SPACEDIM-1]-1;
        int anel_coeff_hi = hi[BL_SPACEDIM-1]+1;

#if (BL_SPACEDIM == 2)
        FORT_SCALEAREA(ax_dat,ARLIM(axlo),ARLIM(axhi), 
                       ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                       anel_coeff[mfi.index()],&anel_coeff_lo,&anel_coeff_hi,lo,hi,&mult);

#elif (BL_SPACEDIM == 3)
        const FArrayBox& zarea = area[2][mfi];
        DEF_CLIMITS(zarea,az_dat,azlo,azhi);
        FORT_SCALEAREA(ax_dat,ARLIM(axlo),ARLIM(axhi), 
                       ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                       az_dat,ARLIM(azlo),ARLIM(azhi), 
                       anel_coeff[mfi.index()],&anel_coeff_lo,&anel_coeff_hi,lo,hi,&mult);

#endif
    }
}

void
MacProj::rescaleArea (int level, MultiFab* area, Real** anel_coeff)
{
    const BoxArray& grids = LevelData[level].boxArray();

    int mult = -1;

    for (MFIter mfi(*area); mfi.isValid(); ++mfi)
    {
        const FArrayBox& xarea  = area[0][mfi];
        const FArrayBox& yarea  = area[1][mfi];

        DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
        DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);

        const int* lo        = grids[mfi.index()].loVect();
        const int* hi        = grids[mfi.index()].hiVect();

        int anel_coeff_lo = lo[BL_SPACEDIM-1]-1;
        int anel_coeff_hi = hi[BL_SPACEDIM-1]+1;

#if (BL_SPACEDIM == 2)
        FORT_SCALEAREA(ax_dat,ARLIM(axlo),ARLIM(axhi), 
                       ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                       anel_coeff[mfi.index()],&anel_coeff_lo,&anel_coeff_hi,lo,hi,&mult);

#elif (BL_SPACEDIM == 3)
        const FArrayBox& zarea = area[2][mfi];
        DEF_CLIMITS(zarea,az_dat,azlo,azhi);
        FORT_SCALEAREA(ax_dat,ARLIM(axlo),ARLIM(axhi), 
                       ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                       az_dat,ARLIM(azlo),ARLIM(azhi), 
                       anel_coeff[mfi.index()],&anel_coeff_lo,&anel_coeff_hi,lo,hi,&mult);

#endif
    }
}
