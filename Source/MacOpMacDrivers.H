
#ifndef _MACOPMACDRIVERS_H_
#define _MACOPMACDRIVERS_H_

#include <AMReX_AmrLevel.H>

#define COVERED_VAL 1.0e40
//====================================================================
// This header file contains declarations for the MacOperator mac projection
// driver routines
//====================================================================

void mlmg_mac_level_solve (amrex::Amr* parent, const amrex::MultiFab* cphi, const amrex::BCRec& phys_bc,
                           const amrex::BCRec& density_math_bc,
                           int level, int Density, amrex::Real mac_tol, amrex::Real mac_abs_tol, amrex::Real rhs_scale,
                           const amrex::MultiFab &S, amrex::MultiFab &Rhs,
                           amrex::MultiFab *u_mac, amrex::MultiFab *mac_phi);

void mlmg_mac_sync_solve (amrex::Amr* parent, const amrex::BCRec& phys_bc,
                          const amrex::BCRec& rho_math_bc,
                          int level, amrex::Real mac_tol, amrex::Real mac_abs_tol, amrex::Real rhs_scale,
                          const amrex::MultiFab& rho, amrex::MultiFab& Rhs,
                          amrex::MultiFab* mac_phi,
			  amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>& Ucorr);

void mlmg_mac_solve (amrex::Amr* parent, const amrex::MultiFab* cphi, const amrex::BCRec& phys_bc,
		     const amrex::BCRec& density_math_bc,
		     int level, amrex::Real mac_tol, amrex::Real mac_abs_tol, amrex::Real rhs_scale,
		     const amrex::MultiFab &rho, amrex::MultiFab &Rhs,
		     amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>& u_mac, amrex::MultiFab *mac_phi,
		     amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>& Ucorr);

#endif
