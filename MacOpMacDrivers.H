#ifndef _MACOPMACDRIVERS_H_

//====================================================================
// This header file contains declarations for the MacOperator mac projection
// driver routines
//====================================================================


// the driver functions for the MAC solves
void mac_level_driver( const MacBndry &mac_bndry,
                       const BoxArray &grids,
                       int use_cg_solve, int level, int Density,
                       const REAL *dx, REAL dt,
                       REAL mac_tol, REAL mac_abs_tol, REAL rhs_scale,
                       MultiFab *area,  MultiFab &volume,
                       MultiFab &S,     MultiFab &Rhs,
                       MultiFab *u_mac, MultiFab *mac_phi );

// a driver function for computing a sync MAC solve
void mac_sync_driver( const MacBndry &mac_bndry,
                      const BoxArray &grids,
                      int use_cg_solve, int level, 
                      const REAL *dx, REAL dt,
                      REAL mac_sync_tol, REAL mac_abs_tol, REAL rhs_scale,
                      MultiFab *area,  MultiFab &volume,
                      MultiFab &Rhs,   MultiFab *rho_half,
                      MultiFab *u_mac, MultiFab *mac_sync_phi );


// apply the mac pressure gradient to a velocity field
// init, means that velocities are initialized here
void mac_vel_update( int init,
                     FARRAYBOX &uvel,
                     FARRAYBOX &vvel,
#if (BL_SPACEDIM == 3 )
                     FARRAYBOX &wvel,
#endif
                     const FARRAYBOX &phi,
                     const FARRAYBOX *rho, int rho_comp,  
                     const BOX &grd, int level, int n,
                     const REAL *dx, REAL scale );

#endif
