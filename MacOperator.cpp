// $Id: MacOperator.cpp,v 1.2 1997-07-17 23:47:15 car Exp $
#include <MacBndry.H>
#include <MacOperator.H>
#include <MacOpMacDrivers.H>
#include <MacOpProjDrivers.H>
#include <MACOPERATOR_F.H>
#include <CGSolver.H>
#include <MultiGrid.H>

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
REAL* fabdat = (fab).dataPtr();
#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const REAL* fabdat = (fab).dataPtr();

// constructor
MacOperator::MacOperator(const BoxArray &ba, const BndryData& mgb, 
                         const REAL *h)
    : ABecLaplacian(ba,mgb,h)
{}



// define the meaning of gradient for the multigrid object
//
void MacOperator::setCoefficients(MultiFab* area, MultiFab& rho, int rho_comp,
			     const REAL* dx)
{
    // should check that all BoxArrays are consistant
    const BoxArray& ba = *gbox[0];
    assert( rho.boxArray() == ba );

    // first set scalar coeficients
    setScalars(0.0,1.0);

    //  NOTE: dont need to set a because alpha is set to zero
    int n_grow = 0;
    MultiFab bxcoef(area[0].boxArray(),area[0].nComp(),n_grow,Fab_allocate);
    MultiFab bycoef(area[1].boxArray(),area[1].nComp(),n_grow,Fab_allocate);
    bxcoef.setVal(0.);
    bycoef.setVal(0.);
#if (BL_SPACEDIM == 3)
    MultiFab bzcoef(area[2].boxArray(),area[2].nComp(),n_grow,Fab_allocate);
    bzcoef.setVal(0.);
#endif

    int ngrd = ba.length();
    int k;
    for (k = 0; k < ngrd; k++) {
	const BOX& grd = ba[k];
        const int* lo = grd.loVect();
        const int* hi = grd.hiVect();
	FArrayBox& bx = bxcoef[k];
	FArrayBox& by = bycoef[k];
	const FArrayBox& ax = area[0][k];
	const FArrayBox& ay = area[1][k];
	const FArrayBox& den = rho[k];

	DEF_LIMITS(bx,bx_dat,bxlo,bxhi);
	DEF_LIMITS(by,by_dat,bylo,byhi);
	DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);
	const int* dlo = den.loVect();
	const int* dhi = den.hiVect();
	const REAL* den_dat = den.dataPtr(rho_comp);

#if (BL_SPACEDIM == 2)
	FORT_MACCOEF(bx_dat,ARLIM(bxlo),ARLIM(bxhi),
                     by_dat,ARLIM(bylo),ARLIM(byhi),
		     ax_dat,ARLIM(axlo),ARLIM(axhi),
                     ay_dat,ARLIM(aylo),ARLIM(ayhi),
		     den_dat,ARLIM(dlo),ARLIM(dhi),lo,hi,dx);
#endif
#if (BL_SPACEDIM == 3)

	FArrayBox& bz = bzcoef[k];
	const FArrayBox& az = area[2][k];

        DEF_CLIMITS(az,az_dat,azlo,azhi);
        DEF_LIMITS(bz,bz_dat,bzlo,bzhi);

	FORT_MACCOEF(bx_dat,ARLIM(bxlo),ARLIM(bxhi),
                     by_dat,ARLIM(bylo),ARLIM(byhi),
                     bz_dat,ARLIM(bzlo),ARLIM(bzhi),
		     ax_dat,ARLIM(axlo),ARLIM(axhi),
                     ay_dat,ARLIM(aylo),ARLIM(ayhi),
                     az_dat,ARLIM(azlo),ARLIM(azhi),
		     den_dat,ARLIM(dlo),ARLIM(dhi),lo,hi,dx);
#endif
    }
  
    bCoefficients(bxcoef,0);
    bCoefficients(bycoef,1);
#if (BL_SPACEDIM == 3)
    bCoefficients(bzcoef,2);
#endif

}


// This function creates the initial rhs for use in the mac multgrid
// solve
//
void MacOperator::defRHS( MultiFab* area, MultiFab& volume, MultiFab& Rhs,
                          MultiFab *vel, REAL scale)
{
      // should check that all BoxArrays are consistant
    const BoxArray& ba = *gbox[0];
    assert( Rhs.boxArray() == ba );

    int ngrd = ba.length();
    int k;
    for (k = 0; k < ngrd; k++) {
	const BOX& grd = ba[k];
        const int* lo = grd.loVect();
        const int* hi = grd.hiVect();
	const FArrayBox& ax = area[0][k];
	const FArrayBox& ay = area[1][k];
	const FArrayBox& vol = volume[k];
	const FArrayBox& ux = vel[0][k];
	const FArrayBox& uy = vel[1][k];
	FArrayBox& rhs = Rhs[k];

	DEF_CLIMITS(ux,ux_dat,uxlo,uxhi);
	DEF_CLIMITS(uy,uy_dat,uylo,uyhi);
	DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);
	DEF_CLIMITS(vol,vol_dat,vlo,vhi);
	DEF_LIMITS(rhs,rhs_dat,rlo,rhi);

#if (BL_SPACEDIM == 2)
	FORT_MACRHS(ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
		    ax_dat,ARLIM(axlo),ARLIM(axhi),
                    ay_dat,ARLIM(aylo),ARLIM(ayhi),
		    vol_dat,ARLIM(vlo),ARLIM(vhi), 
                    rhs_dat,ARLIM(rlo),ARLIM(rhi),
		    lo,hi,&scale);
#endif
#if (BL_SPACEDIM == 3)

	const FArrayBox& az = area[2][k];
        DEF_CLIMITS(az,az_dat,azlo,azhi);

	const FArrayBox& uz = vel[2][k];
	DEF_CLIMITS(uz,uz_dat,uzlo,uzhi);

	FORT_MACRHS(ux_dat,ARLIM(uxlo),ARLIM(uxhi),
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



// apply the mac pressure gradient to a velocity field
// init, means that velocities are initialized here
void mac_vel_update( int init,
                     FARRAYBOX &ux,
                     FARRAYBOX &uy,
#if (BL_SPACEDIM == 3 )
                     FARRAYBOX &uz,
#endif
                     const FARRAYBOX &phi,
                     const FARRAYBOX *rhoptr, int rho_comp,  
                     const BOX &grd, int level, int n,
                     const REAL *dx, REAL scale )
{
    const int* lo = grd.loVect();
    const int* hi = grd.hiVect();
    const FArrayBox& rho = *rhoptr;
    
    DEF_LIMITS(ux,ux_dat,uxlo,uxhi);
    DEF_LIMITS(uy,uy_dat,uylo,uyhi);
    DEF_CLIMITS(phi,phi_dat,p_lo,p_hi);
    const int* rlo = rho.loVect();
    const int* rhi = rho.hiVect();
    const REAL* rho_dat = rho.dataPtr(rho_comp);
    
#if (BL_SPACEDIM == 2)
    FORT_MACUPDATE(&init,
                   ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
                   phi_dat,ARLIM(p_lo),ARLIM(p_hi),
                   rho_dat,ARLIM(rlo),ARLIM(rhi),
                   dx,&scale);
#endif
#if (BL_SPACEDIM == 3)
    
    DEF_LIMITS(uz,uz_dat,uzlo,uzhi);
    
    FORT_MACUPDATE(&init,
                   ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
                   uz_dat,ARLIM(uzlo),ARLIM(uzhi),
                   phi_dat,ARLIM(p_lo),ARLIM(p_hi),
                   rho_dat,ARLIM(rlo),ARLIM(rhi),
                   dx,&scale);
#endif
}


// apply the mac pressure gradient to the divergent mac velocities
// the resultant velocity field is nondivergent
//
void MacOperator::velUpdate(MultiFab* Vel, MultiFab& Phi, const MultiFab& Rho,
                            int rho_comp, const REAL* dx, REAL scale)
{
      // should check that all BoxArrays are consistant
    const BoxArray& ba = *gbox[0];
    assert( Rho.boxArray() == ba );

      // set bndry data in ghost zones
    int apply_lev = 0;
    applyBC(Phi,apply_lev);

    int ngrd = ba.length();
    int k;
    for (k = 0; k < ngrd; k++) {
	const BOX& grd = ba[k];

        mac_vel_update( 0, 
                        Vel[0][k],
                        Vel[1][k],
#if (BL_SPACEDIM == 3 )
                        Vel[2][k],
#endif
                        Phi[k],
                        &(Rho[k]), rho_comp,  
                        grd, 0, k, dx, scale );
    }
}


// multiply by volume*rhs_scale since reflux step (which computed rhs)
// divided by volume.
//
void MacOperator::syncRhs(const MultiFab& Volume, MultiFab& Rhs,
		     REAL rhs_scale, const REAL* dx)
{

    const BoxArray& ba = *gbox[0];

    int ngrd = ba.length();
    int k;
    for (k = 0; k < ngrd; k++) {
	const BOX& grd = ba[k];
        const int* lo = grd.loVect();
        const int* hi = grd.hiVect();
	FArrayBox& rhs = Rhs[k];
	const FArrayBox& vol = Volume[k];

	DEF_CLIMITS(vol,vol_dat,vlo,vhi);
	DEF_LIMITS(rhs,rhs_dat,rlo,rhi);
	FORT_MACSYNCRHS(rhs_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
			vol_dat,ARLIM(vlo),ARLIM(vhi),dx,&rhs_scale);
    }
    Rhs.mult(-1.0,Rhs.nGrow());
}


//========================================================================
//  Driver functions follow
//========================================================================


// a driver function for computing a level MAC solve
void mac_level_driver( const MacBndry &mac_bndry,
                       const BoxArray &grids,
                       int use_cg_solve, int level, int Density,
                       const REAL *dx, REAL dt,
                       REAL mac_tol, REAL mac_abs_tol, REAL rhs_scale,
                       MultiFab *area,  MultiFab &volume,
                       MultiFab &S,     MultiFab &Rhs,
                       MultiFab *u_mac, MultiFab *mac_phi )
{
    MacOperator mac_op(grids,mac_bndry,dx);
    mac_op.setCoefficients(area,S,Density,dx);
    mac_op.defRHS(area,volume,Rhs,u_mac,rhs_scale);
    mac_op.maxOrder(2);
    if (use_cg_solve && mac_op.maxOrder() != 2) {
        cout << "Cant use CGSolver with maxorder > 2 " << endl;
        exit(1);
    }
    
    // construct MultiGrid or CGSolver object and solve system
    if (use_cg_solve) {
        bool use_mg_precond = true;
        CGSolver mac_cg(mac_op,use_mg_precond);
        mac_cg.solve(*mac_phi,Rhs,mac_tol,mac_abs_tol);
    } else {
        MultiGrid mac_mg(mac_op);
        mac_mg.solve(*mac_phi,Rhs,mac_tol,mac_abs_tol);
    }
    
    // velUpdate will set bndry values for mac_phi
    mac_op.velUpdate(u_mac,*mac_phi,S,Density,dx,-dt/2.0);
}


// a driver function for computing a sync MAC solve
void mac_sync_driver( const MacBndry &mac_bndry,
                      const BoxArray &grids,
                      int use_cg_solve, int level, 
                      const REAL *dx, REAL dt,
                      REAL mac_sync_tol, REAL mac_abs_tol, REAL rhs_scale,
                      MultiFab *area,  MultiFab &volume,
                      MultiFab &Rhs,   MultiFab *rho_half,
                      MultiFab *u_mac, MultiFab *mac_sync_phi )
{
    MacOperator mac_op(grids,mac_bndry,dx);
    mac_op.maxOrder(2);
    mac_op.setCoefficients(area,*rho_half, 0, dx);
    mac_op.syncRhs(volume,Rhs,rhs_scale,dx);
    if (use_cg_solve && mac_op.maxOrder() != 2) {
        cout << "Cant use CGSolver with maxorder > 2 " << endl;
        exit(1);
    }
    
    // now construct MultiGrid or CGSolver object to solve system
    if (use_cg_solve) {
        bool use_mg_precond = true;
        CGSolver mac_cg(mac_op,use_mg_precond);
        REAL local_mac_sync_tol = mac_sync_tol/pow(100.0, level);
        mac_cg.solve(*mac_sync_phi,Rhs,local_mac_sync_tol,mac_abs_tol);
    } else {
        MultiGrid mac_mg(mac_op);
        REAL local_mac_sync_tol = mac_sync_tol/pow(100.0, level);
        mac_mg.solve(*mac_sync_phi,Rhs,local_mac_sync_tol,mac_abs_tol);
    }
    
    int mac_op_lev = 0;
    mac_op.applyBC(*mac_sync_phi,mac_op_lev);
}



// make an application specific initial guess
void ProjFirstGuess( MultiFab &U_new, MultiFab &P_new,
                     int level, const BoxArray &grids )
{
}


// scale the variables for a projection solve
void proj_scale_var( MultiFab *rho,         MultiFab *vel,
                     const BoxArray &grids, int level )
{
}

// unscale the variables for a projection solve
void proj_unscale_var( MultiFab *rho,         MultiFab *vel,
                       const BoxArray &grids, int level )
{
}


