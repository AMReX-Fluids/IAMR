
#ifndef _Projection_H_
#define _Projection_H_

#include <BC_TYPES.H>
#include <BCRec.H>
#include <Amr.H>
#include <AmrLevel.H>
#include <RegType.H>
#include <SyncRegister.H>
#include <FluxRegister.H>

class Projection 
{
protected:
  static const int maxlev = 10;
  static const Real BogusValue;
  static const int  LEVEL_PROJ    = 1001;
  static const int  INITIAL_VEL   = 1002;
  static const int  INITIAL_PRESS = 1003;
  static const int  INITIAL_SYNC  = 1004;
  static const int  SYNC_PROJ     = 1005;
  static const int  FILTER_P      = 1006;

public:

    Projection (Amr*   _parent,
                BCRec* _phys_bc, 
                int    _do_sync_proj,
                int    /*_finest_level*/, 
                int    _radius_grow );

    virtual ~Projection ();
    //
    // Set a level in the projector.
    //
    void install_level (int                   level,
                        AmrLevel*             level_data,
                        Array< Array<Real> >* radius);

    void install_anelastic_coefficient (int                        level,
                                        Real                **anel_coeff);
    //
    // Compute the level projection in NavierStokes::advance.
    //
    void level_project (int             level,
                        Real            time,
                        Real            dt,
                        Real            cur_pres_time,
                        Real            prev_pres_time,
                        const Geometry& geom, 
                        MultiFab&       U_old,
                        MultiFab&       U_new,
                        MultiFab&       P_old,
                        MultiFab&       P_new,
                        MultiFab*       rho_half,
                        SyncRegister*   crse_sync_reg,
                        SyncRegister*   fine_sync_reg,
                        int             crse_dt_ratio,
                        int             iteration,
                        int             have_divu);

    // solve DG(correction to P_new) = -D G^perp p^(n-half)
    //  or   DG(correction to P_new) = -D G^perp p^(n-half) - D(U^n /dt)
    void filterP(int level,
                 const Geometry& geom,
                 MultiFab &P_old,
                 MultiFab &P_new,
                 MultiFab &U_old,
                 MultiFab * rho_half,
                 int     ** bc,
                 Real            time,
                 Real            dt,
                 int    have_divu);
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
                      const Geometry& geom,
                      const Real*     dx,
                      Real            dt_crse,
                      int             crse_iteration,
                      int             crse_dt_ratio);
    //
    // A sync project involving the coarse and fine level of a 2 level system.
    //
    void MLsyncProject (int             crse_level,
                        MultiFab&       pres_crse,
                        MultiFab&       vel_crse,
                        MultiFab&       cc_rhs_crse,
                        MultiFab&       pres_fine,
                        MultiFab&       vel_fine,
                        MultiFab&       cc_rhs_fine,
                        MultiFab&       rho_half,
                        MultiFab&       rho_fine,
                        MultiFab*       Vsync,
                        MultiFab&       V_corr,
                        MultiFab&       phi,
                        SyncRegister*   rhs_sync_reg,
                        SyncRegister*   crsr_sync_reg,
                        Real            dt_crse,
                        IntVect&        ratio,
                        int             crse_iteration,
                        int             crse_dt_ratio,
                        const Geometry& crse_geom,
                        bool		pressure_time_is_interval,
                        bool		first_crse_step_after_initial_iters,
                        Real		 cur_crse_pres_time,
                        Real		prev_crse_pres_time,
                        Real		 cur_fine_pres_time,
                        Real		prev_fine_pres_time);
    //
    // This projects an initial vorticity field to define a velocity field.
    //
    void initialVorticityProject (int crse_level);
    //
    // This computes the stream function, given the velocity field
    //  *AND ASSUMING* no-flow boundaries on all physical boundaries.
    //
    void getStreamFunction (PArray<MultiFab>&);
    //
    // The initial velocity projection in post_init.
    // This function ensures that the velocities are nondivergent.
    //
    void initialVelocityProject (int  crse_level,
                                 Real divu_time,
                                 int  have_divu);

    //
    // This function creates an initially hydrostatic pressure field
    //   in the case of nonzero gravity.
    //
    void initialPressureProject (int  crse_level);
    //
    // The velocity projection in post_init, which computes the initial
    // pressure used in the timestepping.
    //
    void initialSyncProject (int       crse_level,
                             MultiFab* sig[],
                             Real      dt, 
                             Real      strt_time,
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
    // Multiply/divide a MultiFab by anel_coeff.
    //
    void AnelCoeffMult (int       level,
                        MultiFab& mf,
                        int       comp);
    void AnelCoeffDiv (int       level,
                       MultiFab& mf,
                       int       comp);

    void put_divu_in_cc_rhs (MultiFab&       rhs,
                             int             level,
                             Real            time);

    //
    // Helper function for NavierStokes::getGradP().
    //
    static void getGradP (FArrayBox& p_fab,
                          FArrayBox& gp,
                          const Box& gpbox_to_fill,
                          const Real* dx);
protected:
    //
    // Read parameters from the input file.
    //
    static void Initialize ();

    static void Finalize ();

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
                 MultiFab&       phi);
    //
    // Convert phi into p^n+1/2.
    //
    void incrPress (int  level,
                    Real dt);
    //
    // This function scales variables at the start of a projection.
    //
    void scaleVar (int             which_call,
                   MultiFab*       sig,
                   int             sig_nghosts,
                   MultiFab*       vel,
                   int             level);
    //
    // This function rescales variables at the end of a projection.
    //
    void rescaleVar (int             which_call,
                     MultiFab*       sig,
                     int             sig_nghosts,
                     MultiFab*       vel,
                     int             level);

    void set_outflow_bcs (int        which_call,
                          MultiFab** phi,
                          MultiFab** Vel_in,
                          MultiFab** Divu_in,
                          MultiFab** Sig_in,
                          int        c_lev,
                          int        f_lev,
                          int        have_divu);

    void set_outflow_bcs_at_level (int          which_call,
                                   int          lev,
                                   int          c_lev,
                                   Box*         state_strip,
                                   Orientation* outFacesAtThisLevel,
                                   int          numOutFlowFaces,
                                   MultiFab**   phi,
                                   MultiFab*    Vel_in,
                                   MultiFab*    Divu_in,
                                   MultiFab*    Sig_in,
                                   int          have_divu,
                                   Real         gravity);

protected:

    void putDown(MultiFab** phi, FArrayBox* phi_fine_strip,
		 int c_lev, int f_lev, const Orientation* outFaces,
		 int numOutFlowFaces, int ncStripWidth);
    //
    // Pointers to amrlevel and amr.
    //
    Amr*             parent;
    PArray<AmrLevel> LevelData;
    //
    // The array of Radii, and number of radii ghost cells.
    //
    int                            radius_grow;
    PArray< Array< Array<Real> > > radius;
    //
    // Anelastic coefficient array.
    //
    Array<Real**> anel_coeff;

    //
    // Boundary objects.
    //
    BCRec*                         phys_bc;

    int do_sync_proj;

    static int  verbose;
    static int  P_code;
    static Real proj_tol;
    static Real sync_tol;
    static Real proj_abs_tol;
    static int  proj_2;
    static int  filter_u;
    static int  rho_wgt_vel_proj;
    static int  do_outflow_bcs;
    static int  make_sync_solvable;
    static Real divu_minus_s_factor;
    static int  add_vort_proj;

    void doNodalProjection(int c_lev, int nlevel, 
			   MultiFab* vel[], MultiFab* phi[], MultiFab* sig[],
			   MultiFab* rhs[], const PArray<MultiFab>& crse_rhs, 
			   Real rel_tol, Real abs_tol,
			   MultiFab* sync_resid_crse=0,
			   MultiFab* sync_resid_fine=0);

  // Mask off fine grids
  // Possible values for msk
  // BogusValue: outside this level
  //          0: covered by fine level
  //          1: valid region not covered fine level
  void mask_grids(MultiFab& msk, const BoxArray& grids, const Geometry& geom,
		  const BoxArray& fineGrids, const IntVect& ref_ratio);

  // Mask off coarse grids
  // Possible values for msk
  // BogusValue: outside physical domain, or we don't care
  //          0: covered by coarse level
  //          1: valid region 
  void mask_grids(MultiFab& msk, const BoxArray& grids, const Geometry& geom);

  // set velocity in ghost cells to zero except for inflow
  void set_boundary_velocity(int c_lev, int nlevel, MultiFab* vel[], bool inflowCorner);

};
#endif
