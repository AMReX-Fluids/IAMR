#ifndef _MacProj_H_
#define _MacProj_H_

//
// $Id: MacProj.H,v 1.5 1997-10-08 20:15:38 car Exp $
//

#include <BC_TYPES.H>
#include <BCRec.H>
#include <Amr.H>
#include <AmrLevel.H>
#include <FluxRegister.H>

class MacProj
{
public:
    // -------------------------------------------------------------
    // Setup functions follow
    // -------------------------------------------------------------
    
    MacProj(Amr * parent,
            int _finest_level,
            BCRec * _phys_bc, 
            int _radius_grow );
    
    virtual ~MacProj();
    void read_params();
    
    void install_level(int level, 
                       AmrLevel * level_data,
                       MultiFab &volume,
                       MultiFab *area,
                       PArray<Real> * _radius );
    
    void setup(int level);
    void cleanup(int level);

    // -------------------------------------------------------------
    // Projection functions follow
    // -------------------------------------------------------------
    
    // the level advance mac projection
    void mac_project(int level, MultiFab * u_mac, MultiFab & S, Real dt,
                     Real prev_time, const MultiFab& divu, int have_divu);
    
    // the sync solve
    void mac_sync_solve(int level, MultiFab * u_mac, Real dt, 
                        MultiFab * rho_half, IntVect& fine_ratio);
    
    // computing the sync tendency
    void mac_sync_compute(int level, MultiFab * u_mac,  
                          MultiFab * Vsync, MultiFab * Ssync,
                          MultiFab * rho_half, 
			  FluxRegister* adv_flux_reg,
                          Array<int> is_conservative,
                          Real prev_time, Real prev_pres_time, Real dt,
                          int NUM_STATE, Real be_cn_theta,
                          const int* increment_sync=NULL);
    
    void mac_sync_compute(int level, MultiFab * u_mac, 
			  MultiFab * Sync, int comp,
                          MultiFab** sync_edges,
			  MultiFab * rho_half,
			  FluxRegister* adv_flux_reg,
			  Array<int> is_conservative, 
			  Real dt);

protected: // ======================== protected member functions

    //------------ test the divergence constraint
    void check_div_cond(int level, MultiFab U_edge[]) const;

    //------------ boundary conditions
    void BuildPhiBC(int level);
    void set_outflow_bcs(int level,
                         MultiFab* mac_phi, Amr* parent, MultiFab* u_mac, 
                         MultiFab& S, const  MultiFab& divu);
    
protected: // ======================== protected data members


    //-------------- pointers to amr,amrlevel
    Amr *parent;
    PArray<AmrLevel> LevelData;
    
    // boundary condition objects
    BCRec * phys_bc;
    Array< Array<BCRec> > phi_bcs;
    
    // MAC sync correction and solution
    PArray<MultiFab> mac_phi_crse;
    PArray<FluxRegister> mac_reg;
    
    // volume and area fractions
    PArray<MultiFab> volume;
    Array<MultiFab*> area;
    
    //-------------- radius data
    int  radius_grow;
    PArray< PArray<Real> > radius;

    //-------------- parameters
    int  finest_level;
    int finest_level_allocated;
    static int verbose;
    static int use_cg_solve;
    static int use_viscosity;
    
    static Real mac_tol;
    static Real mac_abs_tol;
    static Real mac_sync_tol;
    static int do_outflow_bcs;

#if 1
// this code should go away when Marc makes this option
// part of the multigrid code--rbp, 2/13/97
    static int fix_mac_sync_rhs;
#endif
};

#endif
