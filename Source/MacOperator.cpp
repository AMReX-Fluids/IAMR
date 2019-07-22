

#include <AMReX_MacBndry.H>
#include <MacOperator.H>
#include <MacOpMacDrivers.H>
#include <MACOPERATOR_F.H>
#include <AMReX_CGSolver.H>
#include <AMReX_MultiGrid.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();
//
// This is the default maxorder for the linear operator.
//
// To change this ParmParse in a new value for "macop.max_order"
//
static int max_order;

namespace
{
    bool initialized = false;
}

void
MacOperator::Initialize ()
{
    if (initialized) return;
    //
    // This is the default maxorder for the linear operator.
    //
    // To change this ParmParse in a new value for "macop.max_order"
    //
    max_order = 4;

    ParmParse pp("macop");

    pp.query("max_order", max_order);

    amrex::ExecOnFinalize(MacOperator::Finalize);

    initialized = true;
}

void
MacOperator::Finalize ()
{
    initialized = false;
}

MacOperator::MacOperator (Amr*             Parent,
                          const BndryData& mgb,
                          const Real*      h_)
    :
    ABecLaplacian(mgb,h_),
    parent(Parent)
{
    Initialize();
}

MacOperator::~MacOperator () {}

//
// Define the meaning of gradient for the multigrid object.
//

void
MacOperator::setCoefficients (const MultiFab* area,
                              MultiFab&       rho,
                              int             rho_comp,
                              const Real*     dx)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    const BoxArray& ba = gbox[0];
    BL_ASSERT(rho.boxArray() == ba);
    //
    // First set scalar coeficients.
    //
    setScalars(0.0,1.0);
    //
    // Don't need to set a because alpha is set to zero.
    //
    const int n_grow = 0;

    const DistributionMapping& dm = rho.DistributionMap();
    D_TERM(MultiFab bxcoef(area[0].boxArray(),dm,area[0].nComp(),n_grow);,
           MultiFab bycoef(area[1].boxArray(),dm,area[1].nComp(),n_grow);,
           MultiFab bzcoef(area[2].boxArray(),dm,area[2].nComp(),n_grow););
    D_TERM(bxcoef.setVal(0);,
           bycoef.setVal(0);,
           bzcoef.setVal(0););

#ifdef _OPENMP
#pragma omp parallel
#endif	      
    for (MFIter rhomfi(rho,true); rhomfi.isValid(); ++rhomfi)
    {
        BL_ASSERT(ba[rhomfi.index()] == rhomfi.validbox());

        const Box& tilebx    = rhomfi.tilebox();
	// If maccoef were to operate on one dim at a time,
	// then only tilebox would be needed here by switching the MFIter
	// to be on bcoef (which is nodal in one dimension), or perhaps by
	// keeping rhomfi and using enclosedCells/surroudingNodes. 
	// But to keep updating all bcoefs at the same time, must pass validbox
	// and test inside the fortran function 
        const Box& vbx       = rhomfi.validbox();
        const int* lo        = tilebx.loVect();
        const int* hi        = tilebx.hiVect();
	const int* vbxhi     = vbx.hiVect();	
        FArrayBox& bx        = bxcoef[rhomfi];
        FArrayBox& by        = bycoef[rhomfi];
        const FArrayBox& ax  = area[0][rhomfi];
        const FArrayBox& ay  = area[1][rhomfi];
        const FArrayBox& den = rho[rhomfi];

        DEF_LIMITS(bx,bx_dat,bxlo,bxhi);
        DEF_LIMITS(by,by_dat,bylo,byhi);
        DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);

        const int* dlo      = den.loVect();
        const int* dhi      = den.hiVect();
        const Real* den_dat = den.dataPtr(rho_comp);

#if (BL_SPACEDIM == 2)
        maccoef(bx_dat,ARLIM(bxlo),ARLIM(bxhi),
                     by_dat,ARLIM(bylo),ARLIM(byhi),
                     ax_dat,ARLIM(axlo),ARLIM(axhi),
                     ay_dat,ARLIM(aylo),ARLIM(ayhi),
		den_dat,ARLIM(dlo),ARLIM(dhi),lo,hi,vbxhi,dx);
#endif
#if (BL_SPACEDIM == 3)
        FArrayBox& bz       = bzcoef[rhomfi];
        const FArrayBox& az = area[2][rhomfi];

        DEF_CLIMITS(az,az_dat,azlo,azhi);
        DEF_LIMITS(bz,bz_dat,bzlo,bzhi);

        maccoef(bx_dat,ARLIM(bxlo),ARLIM(bxhi),
                     by_dat,ARLIM(bylo),ARLIM(byhi),
                     bz_dat,ARLIM(bzlo),ARLIM(bzhi),
                     ax_dat,ARLIM(axlo),ARLIM(axhi),
                     ay_dat,ARLIM(aylo),ARLIM(ayhi),
                     az_dat,ARLIM(azlo),ARLIM(azhi),
		den_dat,ARLIM(dlo),ARLIM(dhi),lo,hi,vbxhi,dx);
#endif
    }
  
    D_TERM(bCoefficients(bxcoef,0);,
           bCoefficients(bycoef,1);,
           bCoefficients(bzcoef,2););
}

//
// This function creates the initial rhs for use in the mac multgrid solve.
//

void
MacOperator::defRHS (const MultiFab* area,
                     const MultiFab& volume,
                     MultiFab&       Rhs,
                     MultiFab*       vel,
                     Real            scale)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    BL_ASSERT(Rhs.boxArray() == gbox[0]);

#ifdef _OPENMP
#pragma omp parallel
#endif	      
    for (MFIter Rhsmfi(Rhs,true); Rhsmfi.isValid(); ++Rhsmfi)
    {
        const Box& bx        = Rhsmfi.tilebox();
        const int* lo        = bx.loVect();
        const int* hi        = bx.hiVect();
        const FArrayBox& ax  = area[0][Rhsmfi];
        const FArrayBox& ay  = area[1][Rhsmfi];
        const FArrayBox& vol = volume[Rhsmfi];
        const FArrayBox& ux  = vel[0][Rhsmfi];
        const FArrayBox& uy  = vel[1][Rhsmfi];
        FArrayBox& rhs       = Rhs[Rhsmfi];

        DEF_CLIMITS(ux,ux_dat,uxlo,uxhi);
        DEF_CLIMITS(uy,uy_dat,uylo,uyhi);
        DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);
        DEF_CLIMITS(vol,vol_dat,vlo,vhi);
        DEF_LIMITS(rhs,rhs_dat,rlo,rhi);

#if (BL_SPACEDIM == 2)
        macrhs(ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi),
                    ay_dat,ARLIM(aylo),ARLIM(ayhi),
                    vol_dat,ARLIM(vlo),ARLIM(vhi), 
                    rhs_dat,ARLIM(rlo),ARLIM(rhi),
                    lo,hi,&scale);
#endif
#if (BL_SPACEDIM == 3)
        const FArrayBox& az = area[2][Rhsmfi];
        DEF_CLIMITS(az,az_dat,azlo,azhi);

        const FArrayBox& uz = vel[2][Rhsmfi];
        DEF_CLIMITS(uz,uz_dat,uzlo,uzhi);

        macrhs(ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
                    uz_dat,ARLIM(uzlo),ARLIM(uzhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi),
                    ay_dat,ARLIM(aylo),ARLIM(ayhi),
                    az_dat,ARLIM(azlo),ARLIM(azhi),
                    vol_dat,ARLIM(vlo),ARLIM(vhi),
                    rhs_dat,ARLIM(rlo),ARLIM(rhi),
                    lo,hi,&scale);
#endif
    }
    Rhs.mult(-1.0,Rhs.nGrow());
}

//
// Apply the mac pressure gradient to a velocity field.
// init, means that velocities are initialized here.
//

void
mac_vel_update (int              init,
                D_DECL(FArrayBox& ux,
                       FArrayBox& uy,
                       FArrayBox& uz),
                const FArrayBox& phi,
                const FArrayBox& rho,
                int              rho_comp,  
                const Box&       bx,
		const Box&       vbx,
                int              level,
                const Real*      dx,
                Real             scale)
{
    const int* lo        = bx.loVect();
    const int* hi        = bx.hiVect();
    // Need to pass the high coords of validbox into macupdate
    // to test if tilebox is at validbox boundary, in which case
    // we'll need to include an extra point at high end because
    // ux, uy are nodal in one dim (and validbox/phi are cell centered)
    const int* vbx_hi    = vbx.hiVect();
    
    DEF_LIMITS(ux,ux_dat,uxlo,uxhi);
    DEF_LIMITS(uy,uy_dat,uylo,uyhi);
    DEF_CLIMITS(phi,phi_dat,p_lo,p_hi);

    const int* rlo      = rho.loVect();
    const int* rhi      = rho.hiVect();
    const Real* rho_dat = rho.dataPtr(rho_comp);
    
#if (BL_SPACEDIM == 2)
    macupdate(&init,
                   ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
                   phi_dat,ARLIM(p_lo),ARLIM(p_hi),
                   rho_dat,ARLIM(rlo),ARLIM(rhi),
	           lo,hi,vbx_hi,dx,&scale);
#endif
#if (BL_SPACEDIM == 3)
    DEF_LIMITS(uz,uz_dat,uzlo,uzhi);
    
    macupdate(&init,
                   ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
                   uz_dat,ARLIM(uzlo),ARLIM(uzhi),
                   phi_dat,ARLIM(p_lo),ARLIM(p_hi),
                   rho_dat,ARLIM(rlo),ARLIM(rhi),
   	           lo,hi,vbx_hi,dx,&scale);
#endif
}

//
// Apply the mac pressure gradient to the divergent mac velocities.
// The resultant velocity field is nondivergent.
//

void
MacOperator::velUpdate (MultiFab*       Vel,
                        MultiFab&       Phi,
                        const MultiFab& Rho,
                        int             rho_comp,
                        const Real*     dx,
                        Real            scale)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    BL_ASSERT(Rho.boxArray() ==  gbox[0]);
    //
    // Set bndry data in ghost zones.
    //
    int apply_lev = 0;
    applyBC(Phi,0,1,apply_lev);

#ifdef _OPENMP
#pragma omp parallel
#endif	      
    for (MFIter Phimfi(Phi,true); Phimfi.isValid(); ++Phimfi)
    {
        const Box& bx = Phimfi.tilebox();
	const Box& vbx = Phimfi.validbox();
	// If mac_vel_update were to operate on one dim at a time,
	// could pass only one box here, by doing the MFIter on Vel[]
	// (which is nodal in one dimension)
	// and then only tilebox would be needed... Or perhaps could do
	// something smart with surroundingNodes?
	// But to keep updating 
	// all component of Vel at the same time, must pass validbox
	// and test inside the fortran function macupdate
        mac_vel_update(0, 
                       D_DECL(Vel[0][Phimfi],Vel[1][Phimfi],Vel[2][Phimfi]),
                       Phi[Phimfi],
                       Rho[Phimfi], rho_comp,  
                       bx, vbx, 0, dx, scale );
    }
}

//
// Multiply by volume*rhs_scale since reflux step (which computed rhs)
// divided by volume.
//

void
MacOperator::syncRhs (const MultiFab& Volume,
                      MultiFab&       Rhs,
                      Real            rhs_scale,
                      const Real*     dx)
{
  
#ifdef _OPENMP
#pragma omp parallel
#endif	      
  for (MFIter Rhsmfi(Rhs,true); Rhsmfi.isValid(); ++Rhsmfi)
    {
        const Box& grd       = Rhsmfi.tilebox();
        const int* lo        = grd.loVect();
        const int* hi        = grd.hiVect();
        FArrayBox& rhs       = Rhs[Rhsmfi];
        const FArrayBox& vol = Volume[Rhsmfi];

        DEF_CLIMITS(vol,vol_dat,vlo,vhi);
        DEF_LIMITS(rhs,rhs_dat,rlo,rhi);
        macsyncrhs(rhs_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                        vol_dat,ARLIM(vlo),ARLIM(vhi),&rhs_scale);
    }
    Rhs.mult(-1.0,Rhs.nGrow());
}
