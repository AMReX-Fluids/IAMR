//BL_COPYRIGHT_NOTICE

//
// $Id: Diffusion.H,v 1.10 1999-02-18 18:22:14 sstanley Exp $
//

#ifndef _Diffusion_H_
#define _Diffusion_H_

//
// Comment out this line to use diffusion class outside
// the context of NavierStokes and classes derived from it.
//
#define USE_NAVIERSTOKES 1

#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LevelBld.H>
#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <ABecLaplacian.H>
#include <ViscBndry.H>

//
// Include files for tensor solve.
//
#include <ViscBndryTensor.H>
#include <DivVis.H>
#include <LO_BCTYPES.H>
#include <MCMultiGrid.H>
#include <MCCGSolver.H>

#ifndef _NavierStokes_H_
enum StateType {State_Type=0, Press_Type}; 
#  if (BL_SPACEDIM == 2)
enum StateNames  { Xvel=0, Yvel, Density};
#  else
enum StateNames  { Xvel=0, Yvel, Zvel, Density};
#  endif
#endif

class Diffusion 
{
public:

    enum SolveMode {PREDICTOR, CORRECTOR, ONEPASS};

    Diffusion ();

    Diffusion (Amr*               Parent,
               AmrLevel*          Caller,
               Diffusion*         coarser,
               int                num_state,
               FluxRegister*      Viscflux_reg,
               MultiFab&          Volume,
               MultiFab*          Area,
               const Array<int>&  _is_diffusive,
               const Array<Real>& _visc_coef);

    virtual ~Diffusion ();

    void echo_settings () const;

    FluxRegister* viscFluxReg () { return viscflux_reg; }

    REAL get_scaled_abs_tol (int                    sigma,
                             const MultiFab*        rhs,
                             REAL                   a,
                             REAL                   b,
                             const MultiFab*        alpha,
                             const MultiFab* const* beta,
                             REAL                   reduction) const;

    REAL get_scaled_abs_tol (int                    sigma,
                             const MultiFab*        rhs,
                             REAL                   a,
                             REAL                   b,
                             const MultiFab*        alpha,
                             const MultiFab* const* betan,
                             const MultiFab* const* betanp1,
                             REAL                   reduction) const;

    void diffuse_scalar (Real                   dt,
			 int                    sigma,
			 Real                   be_cn_theta,
			 const MultiFab*        rho_half,
			 int                    rho_flag,
			 MultiFab* const*       flux,
			 int                    dataComp = 0,
			 MultiFab*              delta_rhs = 0, 
			 const MultiFab*        alpha = 0, 
			 const MultiFab* const* betan = 0, 
			 const MultiFab* const* betanp1 = 0,
			 const SolveMode&       solve_mode = ONEPASS);
    
    void diffuse_velocity (Real                   dt,
                           Real                   be_cn_theta,
                           const MultiFab*        rho_half,
                           int                    rho_flag,
                           MultiFab*              delta_rhs = 0, 
                           const MultiFab* const* betan = 0, 
                           const MultiFab* const* betanp1 = 0);

    void diffuse_velocity_constant_mu (Real            dt,
                                       Real            be_cn_theta,
                                       const MultiFab* rho_half,
                                       MultiFab*       delta_rhs);

    void diffuse_tensor_velocity (Real                   dt,
                                  Real                   be_cn_theta,
                                  const MultiFab*        rho_half,
                                  MultiFab*              delta_rhs, 
                                  const MultiFab* const* betan, 
                                  const MultiFab* const* betanp1);

    void diffuse_Vsync (MultiFab*              Vsync,
                        Real                   dt,
                        Real                   be_cn_theta,
                        const MultiFab*        rho_half,
                        int                    rho_flag,
                        const MultiFab* const* beta = 0);

    void diffuse_Vsync_constant_mu (MultiFab*       Vsync,
                                    Real            dt,
                                    Real            be_cn_theta,
                                    const MultiFab* rho_half,
                                    int             rho_flag);

    void diffuse_tensor_Vsync (MultiFab*              Vsync,
                               Real                   dt,
                               Real                   be_cn_theta,
                               const MultiFab*        rho_half,
                               int                    rho_flag,
                               const MultiFab* const* beta);

    void diffuse_Ssync (MultiFab*              Ssync,
                        int                    sigma,
                        Real                   dt,
                        Real                   be_cn_theta,
                        const MultiFab*        rho_half,
                        int                    rho_flag,
			MultiFab* const*       flux,
			int                    dataComp = 0,
                        const MultiFab* const* beta = 0,
			const MultiFab*        alpha = 0);

    ABecLaplacian* getViscOp (int                    src_comp,
                              const Real             a,
                              const Real             b, 
                              const Real             time,
                              ViscBndry&             visc_bndry,
                              const MultiFab*        rho_half,
                              int                    rho_flag,
                              Real*                  rhsscale = 0,
			      int                    dataComp = 0,
                              const MultiFab* const* beta = 0,
                              const MultiFab*        alpha_in = 0);

    ABecLaplacian* getViscOp (int                    src_comp,
                              const Real             a,
                              const Real             b,
                              const MultiFab*        rho_half,
                              int                    rho_flag,
                              Real*                  rhsscale = 0,
			      int                    dataComp = 0,
                              const MultiFab* const* beta = 0,
                              const MultiFab*        alpha_in = 0);

    DivVis* getTensorOp (Real                   a,
                         Real                   b, 
                         Real                   time, 
                         ViscBndryTensor&       visc_bndry,
                         const MultiFab*        rho_half,
			 int                    dataComp,
                         const MultiFab* const* beta);

    void getTensorBndryData(
        ViscBndryTensor& bndry, 
        Real time);

    DivVis* getTensorOp (Real                   a,
                         Real                   b, 
                         const MultiFab*        rho_half,
			 int                    dataComp,
                         const MultiFab* const* beta);

    void getViscTerms (MultiFab&              visc_terms,
                       int                    src_comp,
                       int                    comp,
                       Real                   time,
                       int                    rho_flag,
		       int                    dataComp = 0,
                       const MultiFab* const* beta = 0);

    void getTensorViscTerms (MultiFab&              visc_terms, 
                             Real                   time,
			     int                    dataComp,
                             const MultiFab* const* beta);

    void getBndryData (ViscBndry& bndry,
                       int        src_comp,
                       int        num_comp,
                       Real       time,
                       int        rho_flag);

    void FillBoundary (BndryRegister& bdry,
                       int            src_comp,
                       int            dest_comp,
                       int            num_comp,
                       Real           time,
                       int            rho_flag);

    int useDvConstantMu () { return use_dv_constant_mu_def; }

    void checkBetas (const MultiFab* const* beta1,
                     const MultiFab* const* beta2,
                     int&                   allthere,
                     int&                   allnull) const;

    void checkBeta (const MultiFab* const* beta,
                    int&                   allthere,
                    int&                   allnull) const;

    void checkBeta (const MultiFab* const* beta,
                    int&                   allthere) const;

    void allocFluxBoxesLevel (MultiFab**& fluxbox, 
                              int         nghost = 0,
                              int         nvar = 1);

    void removeFluxBoxesLevel (MultiFab**& fluxbox);

#ifdef USE_NAVIERSTOKES
    void compute_divmusi (Real                   time,
			  const MultiFab* const* beta,
			  MultiFab&              divmusi);

    void compute_divmusi (Real      time,
			  Real      mu,
			  MultiFab& divmusi);
#endif

    inline int maxOrder() const;
    inline int tensorMaxOrder() const;

private:
    //
    // The data.
    //
    Amr*            parent;
    AmrLevel*       caller;
    Diffusion*      coarser;
    Diffusion*      finer;
    const BoxArray& grids;
    int             level;
    int             NUM_STATE;
    IntVect         crse_ratio;
    FluxRegister*   viscflux_reg;
    //
    // Volume and area fractions.
    //
    MultiFab& volume;
    MultiFab* area;
    //
    // Static data.
    //
    static int         first;
    static int         do_reflux;
    static int         use_cg_solve;
    static int         use_tensor_cg_solve;
    static int         max_order;
    static int         tensor_max_order;
    static bool        use_mg_precond_flag;
    static int         use_dv_constant_mu_def;
    static int         scale_abec;
    static int         est_visc_mag;
    static int         Lphi_in_abs_tol;
    static Array<int>  is_diffusive;    // Does variable diffuse?
    static Array<Real> visc_coef;       // Const coef viscosity terms
    static Real        visc_tol;
    static Real        visc_abs_tol;
    static int         verbose;
    static Array<REAL> typical_vals;
};

inline
REAL
Diffusion::get_scaled_abs_tol (int             sigma,
                               const MultiFab* rhs,
                               REAL            a,
                               REAL            b,
                               const MultiFab* alpha,
                               const MultiFab* const * beta,
                               REAL            reduction) const
{
    return get_scaled_abs_tol(sigma, rhs, a, b, alpha, beta, 0, reduction);
}

inline
int
Diffusion::maxOrder() const
{
    return max_order;
}

inline
int
Diffusion::tensorMaxOrder() const
{
    return tensor_max_order;
}

#endif

