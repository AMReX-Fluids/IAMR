
#ifndef _MacProj_H_
#define _MacProj_H_

#include <BC_TYPES.H>
#include <BCRec.H>
#include <Amr.H>
#include <AmrLevel.H>
#include <FluxRegister.H>
#include <Godunov.H>

class MacProj
{
public:
    
    MacProj (Amr*   parent,
             int    _finest_level,
             BCRec* _phys_bc, 
             int    /*not used*/ );
    
    virtual ~MacProj ();

    
    void install_level (int       level, 
                        AmrLevel* level_data);

    void install_anelastic_coefficient (int                    level, 
                                        Real**            anel_coeff);
    
    void setup (int level);
    void cleanup (int level);
    //
    // The level advance mac projection.
    //
    void mac_project (int             level,
                      MultiFab*       u_mac,
                      MultiFab&       S,
                      Real            dt,
                      Real            prev_time,
                      const MultiFab& divu,
                      int             have_divu,
                      bool            increment_vel_register = true);

    void test_umac_periodic (int level, MultiFab* u_mac);
    //
    // The sync solve.
    //
    void mac_sync_solve (int       level,
                         Real      dt, 
                         MultiFab* rho_half,
                         IntVect&  fine_ratio);
    //
    // Computing the sync tendency.
    //
    void mac_sync_compute (int                   level,
                           MultiFab*             u_mac, 
                           MultiFab*             Vsync,
                           MultiFab*             Ssync,
                           MultiFab*             rho_half, 
                           FluxRegister*         adv_flux_reg,
                           Array<AdvectionForm>& advectionType,
                           Real                  prev_time,
                           Real                  prev_pres_time,
                           Real                  dt,
                           int                   NUM_STATE,
                           Real                  be_cn_theta,
                           bool                  modify_reflux_normal_vel,
                           int		         do_mom_diff,
                           const Array<int>&     increment_sync = Array<int>());

    void mac_sync_compute (int                   level,
                           MultiFab*             Sync,
                           int                   comp,
                           int                   s_ind,
                           const MultiFab*       const* sync_edges,
			   int                   eComp,
                           MultiFab*             rho_half,
                           FluxRegister*         adv_flux_reg,
                           Array<AdvectionForm>& advectionType,
                           bool		         modify_reflux_normal_vel,
                           Real                  dt);

    void   scaleArea (int level, MultiFab* area, Real** anel_coeff);
    void rescaleArea (int level, MultiFab* area, Real** anel_coeff);

    static void Initialize ();
    static void Finalize ();

protected:
    //
    // Test the divergence constraint.
    //
    void check_div_cond (int      level,
                         MultiFab U_edge[]) const;
    //
    // Boundary conditions.
    //
    void BuildPhiBC (int level);

    void set_outflow_bcs (int level,
                          MultiFab* mac_phi,
                          const MultiFab* u_mac, 
                          const MultiFab& S,
                          const  MultiFab& divu);
    //
    // Pointers to amr,amrlevel.
    //
    Amr*             parent;
    PArray<AmrLevel> LevelData;
    //
    // Boundary condition objects.
    //
    BCRec*                phys_bc;
    Array< Array<BCRec> > phi_bcs;
    //
    // MAC sync correction and solution.
    //
    PArray<MultiFab>     mac_phi_crse;
    PArray<FluxRegister> mac_reg;
    //
    // Anelastic coefficient array.
    //
    Array<Real**> anel_coeff;
    //
    // Parameters.
    //
    int        finest_level;
    int        finest_level_allocated;

    static int  verbose;
    static bool  use_cg_solve;
    static int  do_outflow_bcs;
    static Real mac_tol;
    static Real mac_abs_tol;
    static Real mac_sync_tol;
    static int  check_umac_periodicity;
    static int fix_mac_sync_rhs;
};

#endif
