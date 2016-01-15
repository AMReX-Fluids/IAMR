
#ifndef _Diffusion_H_
#define _Diffusion_H_

#include <winstd.H>

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
#include <FluxBoxes.H>

//
// Include files for tensor solve.
//
#include <ViscBndryTensor.H>
#include <DivVis.H>
#include <LO_BCTYPES.H>
#include <MCMultiGrid.H>
#include <MCCGSolver.H>

//
// Useful enumeration of the different forms of the diffusion terms
//
enum DiffusionForm { RhoInverse_Laplacian_S, Laplacian_SoverRho, Laplacian_S };

class NavierStokesBase;

class Diffusion 
{
public:

    enum SolveMode {PREDICTOR, CORRECTOR, ONEPASS};

    Diffusion ();

    Diffusion (Amr*               Parent,
               NavierStokesBase*  Caller,
               Diffusion*         coarser,
               int                num_state,
               FluxRegister*      Viscflux_reg,
               const Array<int>&  _is_diffusive,
               const Array<Real>& _visc_coef);

    ~Diffusion ();

    void echo_settings () const;

    FluxRegister* viscFluxReg ();

    Real get_scaled_abs_tol (const MultiFab& rhs,
                             Real            reduction) const;

    void diffuse_scalar (Real                   dt,
			 int                    sigma,
			 Real                   be_cn_theta,
			 const MultiFab&        rho_half,
			 int                    rho_flag,
			 MultiFab* const*       fluxn,
			 MultiFab* const*       fluxnp1,
			 int                    fluxComp,
			 MultiFab*              delta_rhs, 
                         int                    rhsComp,
			 const MultiFab*        alpha, 
                         int                    alphaComp,
			 const MultiFab* const* betan, 
			 const MultiFab* const* betanp1,
                         int                    betaComp,
			 const SolveMode&       solve_mode = ONEPASS,
                         bool                   add_old_time_divFlux = true);
    
    void diffuse_velocity (Real                   dt,
                           Real                   be_cn_theta,
                           const MultiFab&        rho_half,
                           int                    rho_flag,
                           MultiFab*              delta_rhs, 
                           int                    rhsComp,
                           const MultiFab* const* betan, 
                           const MultiFab* const* betanp1,
                           int                    betaComp);

    void diffuse_velocity (Real                   dt,
                           Real                   be_cn_theta,
                           const MultiFab&        rho_half,
                           int                    rho_flag,
                           MultiFab*              delta_rhs, 
                           const MultiFab* const* betan, 
                           const MultiFab* const* betanp1);

    void diffuse_velocity_constant_mu (Real            dt,
                                       Real            be_cn_theta,
                                       const MultiFab& rho_half,
                                       MultiFab*       delta_rhs,
                                       int             rhsComp);

    void diffuse_tensor_velocity (Real                   dt,
                                  Real                   be_cn_theta,
                                  const MultiFab&        rho_half,
                                  int                    rho_flag,
                                  MultiFab*              delta_rhs, 
                                  int                    rhsComp,
                                  const MultiFab* const* betan, 
                                  const MultiFab* const* betanp1,
                                  int                    betaComp);

    void diffuse_Vsync (MultiFab&              Vsync,
                        Real                   dt,
                        Real                   be_cn_theta,
                        const MultiFab&        rho_half,
                        int                    rho_flag,
                        const MultiFab* const* beta,
                        int                    betaComp = 0);

    void diffuse_Vsync_constant_mu (MultiFab&       Vsync,
                                    Real            dt,
                                    Real            be_cn_theta,
                                    const MultiFab& rho_half,
                                    int             rho_flag);

    void diffuse_tensor_Vsync (MultiFab&              Vsync,
                               Real                   dt,
                               Real                   be_cn_theta,
                               const MultiFab&        rho_half,
                               int                    rho_flag,
                               const MultiFab* const* beta,
                               int                    betaComp);

    void diffuse_Ssync (MultiFab&              Ssync,
                        int                    sigma,
                        Real                   dt,
                        Real                   be_cn_theta,
                        const MultiFab&        rho_half,
                        int                    rho_flag,
			MultiFab* const*       flux,
			int                    fluxComp,
                        const MultiFab* const* beta,
                        int                    betaComp,
			const MultiFab*        alpha,
                        int                    alphaComp);

    ABecLaplacian* getViscOp (int                    src_comp,
                              Real                   a,
                              Real                   b, 
                              Real                   time,
                              ViscBndry&             visc_bndry,
                              const MultiFab&        rho_half,
                              int                    rho_flag,
                              Real*                  rhsscale,
                              const MultiFab* const* beta,
			      int                    betaComp,
                              const MultiFab*        alpha_in,
			      int                    alphaComp,
                              bool		     bndry_already_filled = false);

    ABecLaplacian* getViscOp (int                    src_comp,
                              Real                   a,
                              Real                   b,
                              const MultiFab&        rho_half,
                              int                    rho_flag,
                              Real*                  rhsscale,
                              const MultiFab* const* beta,
			      int                    betaComp,
                              const MultiFab*        alpha_in,
			      int                    alphaComp);

    void setAlpha (ABecLaplacian*  visc_op,
                   int             src_comp,
                   Real            a,
                   Real            b, 
                   Real            time,
                   const MultiFab& rho_half,
                   int                    rho_flag,
                   Real*                  rhsscale,
                   int                    dataComp,
                   const MultiFab* alpha_in);

    void setBeta (ABecLaplacian*         visc_op,
                  const MultiFab* const* beta,
                  int                    betaComp);

    DivVis* getTensorOp (Real                   a,
                         Real                   b, 
                         Real                   time, 
                         ViscBndryTensor&       visc_bndry,
                         const MultiFab&        rho_half,
                         const MultiFab* const* beta,
			 int                    betaComp);

    void getTensorBndryData (ViscBndryTensor& bndry, Real time);

    DivVis* getTensorOp (Real                   a,
                         Real                   b, 
                         const MultiFab&        rho_half,
                         const MultiFab* const* beta,
			 int                    betaComp);

    void getViscTerms (MultiFab&              visc_terms,
                       int                    src_comp,
                       int                    comp,
                       Real                   time,
                       int                    rho_flag,
                       const MultiFab* const* beta,
		       int                    betaComp);

    void getTensorViscTerms (MultiFab&              visc_terms, 
                             Real                   time,
                             const MultiFab* const* beta,
			     int                    dataComp);

    void getBndryData (ViscBndry&         bndry,
                       int                state_ind,
                       int                num_comp,
                       Real               time,
                       int                rho_flag);

    void getBndryData (ViscBndry&          bndry,
                       int                 state_ind,
                       int                 num_comp,
                       AmrLevel::TimeLevel the_time,
                       int                 rho_flag);

    void getBndryDataGivenS (ViscBndry&         bndry,
                             MultiFab&          S,
                             MultiFab&          S_crse,
                             int                state_ind,
                             int                src_comp,
                             int                num_comp);

    void FillBoundary (BndryRegister&     bdry,
                       int                src_comp,
                       int                dest_comp,
                       int                num_comp,
                       Real               time,
                       int                rho_flag);

    void checkBetas (const MultiFab* const* beta1,
                     const MultiFab* const* beta2,
                     int&                   allthere,
                     int&                   allnull) const;

    void checkBeta (const MultiFab* const* beta,
                    int&                   allthere,
                    int&                   allnull) const;

    void checkBeta (const MultiFab* const* beta,
                    int&                   allthere) const;

    void compute_divmusi (Real                   time,
			  const MultiFab* const* beta,
			  MultiFab&              divmusi);

    void compute_divmusi (Real      time,
			  Real      mu,
			  MultiFab& divmusi);

    int maxOrder () const;
    int tensorMaxOrder () const;

    static int set_rho_flag (const DiffusionForm compDiffusionType);

    static bool are_any (const Array<DiffusionForm>& diffusionType,
                         const DiffusionForm         testForm,
                         const int                   sComp,
                         const int                   nComp);

    static int how_many (const Array<DiffusionForm>& diffusionType,
                         const DiffusionForm         testForm,
                         const int                   sComp,
                         const int                   nComp);

protected:

    void getTensorOp_doit (DivVis*                tensor_op,
                           Real                   a,
                           Real                   b,
                           const MultiFab&        rho,
                           const MultiFab* const* beta,
                           int                    betaComp);

    static void Finalize ();
    //
    // Data Required by Derived Classes
    //
    Amr*              parent;
    NavierStokesBase* navier_stokes;
    const BoxArray&   grids;
    const int         level;
    //
    // Static data.
    //
    static bool        use_mg_precond_flag;
    static int         use_cg_solve;
    static int         scale_abec;
    static Array<int>  is_diffusive;    // Does variable diffuse?
    static Array<Real> visc_coef;       // Const coef viscosity terms
    static int         verbose;
    static Real        visc_tol;

private:
    //
    // The data.
    //
    Diffusion*    coarser;
    Diffusion*    finer;
    int           NUM_STATE;
    IntVect       crse_ratio;
    FluxRegister* viscflux_reg;
    //
    // Static data.
    //
    static int         do_reflux;
    static int         use_tensor_cg_solve;
    static int         max_order;
    static int         tensor_max_order;
};

#endif

