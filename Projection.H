//BL_COPYRIGHT_NOTICE

#ifndef _Projection_H_
#define _Projection_H_

//
// $Id: Projection.H,v 1.16 1998-11-05 23:30:18 lijewski Exp $
//

#include <BC_TYPES.H>
#include <BCRec.H>
#include <Amr.H>
#include <AmrLevel.H>
#include <RegType.H>
#include <SyncRegister.H>
#include <FluxRegister.H>

class holy_grail_amr_projector;
class inviscid_fluid_boundary_class;

class Projection 
{
public:

    Projection (Amr*   _parent,
                BCRec* _phys_bc, 
                int    _do_sync_proj,
                int    _finest_level, 
                int    _radius_grow );

    virtual ~Projection ();
    //
    // Set a level in the projector.
    //
    void install_level (int           level,
                        AmrLevel*     level_data,
                        PArray<Real>* _radius );
    //
    // Reset finest_level.
    //
    void setFinestLevel (int new_finest) { finest_level = new_finest; }
    //
    // Compute the level projection in NavierStokes::advance.
    //
    void level_project (int             level,
                        Real            dt,
                        Real            cur_pres_time,
                        Real            old_state_time,
                        Real            cur_state_time,
                        const Geometry& geom, 
                        MultiFab&       U_old,
                        MultiFab&       U_new,
                        MultiFab&       P_old,
                        MultiFab&       P_new,
                        MultiFab*       rho_half,
                        MultiFab&       dsdt, 
                        SyncRegister*   crse_sync_reg,
                        SyncRegister*   fine_sync_reg,
                        int             crse_dt_ratio,
                        int**           bc,
                        int             iteration,
                        Real            divu_minus_s_factor,
                        MultiFab&       divuold,
                        int             have_divu);
    //
    // Attempt at creating better pressure fill.
    //
    void harmonic_project (int             level,
                           Real            dt,
                           Real            cur_pres_time,
                           const Geometry& geom,
                           MultiFab&       P_old);
    //
    // A sync project involving only the coarse level of a 2 level system.
    //
    void syncProject (int             crse_level,
                      MultiFab&       pres,
                      MultiFab&       vel,
                      MultiFab*       rho_half,
                      MultiFab*       Vsync,
                      MultiFab&       phi,
                      SyncRegister*   rhs_sync_reg,
                      SyncRegister*   crsr_sync_reg, 
                      const BoxArray& sync_boxes,
                      int**           bc,
                      const Geometry& geom,
                      const Real*     dx,
                      Real            dt_crse,
                      int             crse_dt_ratio);
    //
    // A sync project involving the coarse and fine level of a 2 level system.
    //
    void MLsyncProject (int             crse_level,
                        MultiFab&       pres_crse,
                        MultiFab&       vel_crse,
                        MultiFab&       pres_fine,
                        MultiFab&       vel_fine,
                        MultiFab&       rho_half,
                        MultiFab&       rho_fine,
                        MultiFab*       Vsync,
                        MultiFab&       V_corr,
                        MultiFab&       phi,
                        SyncRegister*   rhs_sync_reg,
                        SyncRegister*   crsr_sync_reg,
                        int**           bc,
                        const Real*     dx,
                        Real            dt_crse,
                        IntVect&        ratio,
                        int             crse_dt_ratio,
                        const Geometry& fine_geom,
                        const Geometry& crse_geom);
    //
    // This projects an initial vorticity field to define a velocity field.
    //
    void initialVorticityProject (int crse_level);
    //
    // The initial velocity projection in post_init.
    // This function ensures that the velocities are nondivergent.
    //
    void initialVelocityProject (int  crse_level,
                                 Real divu_time,
                                 int  have_divu);
    //
    // The velocity projection in post_init, which computes the initial
    // pressure used in the timestepping.
    //
    void initialSyncProject (int       crse_level,
                             MultiFab* sig[],
                             Real      dt, 
                             Real      strt_time,
                             Real      dt_init,
                             int       have_divu);
    //
    // Enforce periodicity on a multifab.
    //
    void EnforcePeriodicity (MultiFab&       psi,
                             int             nvar,
                             const BoxArray& grids,
                             const Geometry& geom);
    //
    // Multiply/divide a MultiFab by radius for r-z coordinates.
    // These should really be protected.
    //
    void radMult (int       level,
                  MultiFab& mf,
                  int       comp);

    void radDiv (int       level,
                 MultiFab& mf,
                 int       comp);
    //
    // Filter the velocities and pressures after a projection.
    //
    void filterUandP (int             level,
                      MultiFab&       P_new, 
                      MultiFab&       U_new,
                      MultiFab*       rho_half,
                      const BoxArray& grids,
                      const Real*     dx,
                      const Real      dt);
    //
    // Compute the node-centered divergence of U.
    //
    void computeDV (MultiFab&       DV,
                    const MultiFab& U,
                    int             src_comp,
                    const Real*     dx,
                    int             is_rz);
    //
    // Put divu in the rhs of the projector.
    //
    void put_divu_in_node_rhs (MultiFab&       rhs,
                               int             level,
                               const int&      nghost,
                               const BoxArray& P_grids,
                               Real            time,
                               int             user_rz = -1);

    void put_divu_in_cc_rhs (MultiFab&       rhs,
                             int             level,
                             const BoxArray& grids,
                             Real            time);

protected:
    //
    // Read parameters from the input file.
    //
    static void read_params ();
    //
    // Set up the sync projection.
    //
    void bldSyncProject ();
    //
    // Fill patch to get the divU.
    //
    void getDivCond (int        level,
                     FArrayBox& fab,
                     int        ngrow,
                     Real       time);

    MultiFab* getDivCond (int  level,
                          int  ngrow,
                          Real time);
    //
    // Convert U from an Accleration like quantity to a velocity
    // Unew = Uold + alpha*Unew
    //
    void UnConvertUnew (MultiFab&       Uold,
                        Real            alpha,
                        MultiFab&       Unew, 
                        const BoxArray& grids );
    //
    // Convert U from an Accleration like quantity to a velocity
    // Unew = Uold + alpha*Unew
    //
    void UnConvertUnew (FArrayBox& Uold,
                        Real       alpha,
                        FArrayBox& Unew,
                        const Box& grd);
    //
    // Convert U to an Accleration like quantity
    // Unew = (Unew - Uold)/alpha
    //
    void ConvertUnew (MultiFab&       Unew,
                      MultiFab&       Uold,
                      Real            alpha,
                      const BoxArray& grids );
    //
    // Convert U to an Accleration like quantity.
    // Unew = (Unew - Uold)/alpha
    //
    void ConvertUnew (FArrayBox& Unew,
                      FArrayBox& Uold,
                      Real       alpha,
                      const Box& grd);
    //
    // Update a quantity U using the formula
    // Unew = Unew + alpha*Uold
    //
    void UpdateArg1 (MultiFab&       Unew,
                     Real            alpha,
                     MultiFab&       Uold,
                     int             nvar,
                     const BoxArray& grids,
                     int ngrow);
    //
    // Update a quantity U using the formula
    // Unew = Unew + alpha*Uold
    //
    void UpdateArg1 (FArrayBox& Unew,
                     Real       alpha,
                     FArrayBox& Uold,
                     int        nvar,
                     const Box& grd,
                     int        ngrow);
    //
    // Add phi to P.
    //
    void AddPhi (MultiFab&       p,
                 MultiFab&       phi,
                 const BoxArray& grids);
    //
    // Convert phi into p^n+1/2.
    //
    void incrPress (int  level,
                    Real dt);
    //
    // This function scales variables at the start of a projection.
    //
    void scaleVar (MultiFab*       sig,
                   int             sig_nghosts,
                   MultiFab*       vel,
                   const BoxArray& grids,
                   int             level);
    //
    // This function rescales variables at the end of a projection.
    //
    void rescaleVar (MultiFab*       sig,
                     int             sig_nghosts,
                     MultiFab*       vel,
                     const BoxArray& grids,
                     int             level);

    void setUpBcs ();

    void set_level_projector_outflow_bcs (int       level,
                                          Real      cur_state_time,
                                          MultiFab& P_new);
    
    void set_initial_projection_outflow_bcs (MultiFab** vel,
                                             MultiFab** sig,
                                             MultiFab** phi,
                                             int        c_lev,
                                             Real       cur_divu_time);
    
    void set_initial_syncproject_outflow_bcs (MultiFab** phi,
                                              int        c_lev,
                                              Real       start_time,
                                              Real       dt);
protected:
    //
    // Pointers to amrlevel and amr.
    //
    Amr*             parent;
    PArray<AmrLevel> LevelData;
    //
    // The array of Radii, and number of radii ghost cells.
    //
    int                    radius_grow;
    PArray< PArray<Real> > radius;
    //
    // The aliasLib level projector.
    //
    holy_grail_amr_projector *sync_proj;
    //
    // Boundary objects.
    //
    BCRec*                         phys_bc;
    inviscid_fluid_boundary_class* projector_bndry;

    int finest_level;
    int do_sync_proj;

    static int  verbose;
    static int  P_code;
    static Real proj_tol;
    static Real sync_tol;
    static Real proj_abs_tol;
    static Real filter_factor;
    static int  filter_u;
    static int  rho_wgt_vel_proj;
    static int  do_outflow_bcs;
    static int  make_sync_solvable;
};
#endif
