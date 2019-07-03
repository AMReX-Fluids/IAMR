
#ifndef _MACOPMACDRIVERS_H_
#define _MACOPMACDRIVERS_H_

#include <MacOperator.H>

//====================================================================
// This header file contains declarations for the MacOperator mac projection
// driver routines
//====================================================================


// Apply the mac pressure gradient to a velocity field 
// "init" means that velocities are initialized here
void mac_vel_update( int init,
                     amrex::FArrayBox &uvel,
                     amrex::FArrayBox &vvel,
#if (BL_SPACEDIM == 3 )
                     amrex::FArrayBox &wvel,
#endif
                     const amrex::FArrayBox &phi,
                     const amrex::FArrayBox &rho, int rho_comp,  
                     const amrex::Box &bx, 
		     const amrex::Box &vbx, int level,
                     const amrex::Real *dx, amrex::Real scale );

#endif


void mlmg_mac_level_solve (amrex::Amr* parent, const amrex::MultiFab* cphi, const amrex::BCRec& phys_bc,
                           int level, int Density, amrex::Real mac_tol, amrex::Real mac_abs_tol, amrex::Real rhs_scale,
                           const amrex::MultiFab *area, const amrex::MultiFab &volume,
                           const amrex::MultiFab &S, amrex::MultiFab &Rhs,
                           amrex::MultiFab *u_mac, amrex::MultiFab *mac_phi, int verbose);

void mlmg_mac_sync_solve (amrex::Amr* parent, const amrex::BCRec& phys_bc,
                          int level, amrex::Real mac_tol, amrex::Real mac_abs_tol, amrex::Real rhs_scale,
                          const amrex::MultiFab* area, const amrex::MultiFab& volume,
                          const amrex::MultiFab& rho, amrex::MultiFab& Rhs,
                          amrex::MultiFab* mac_phi, int verbose);
