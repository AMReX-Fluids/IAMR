
#ifndef _Diffusion_H_
#define _Diffusion_H_


#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_LevelBld.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_ErrorList.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_ABecLaplacian.H>
#include <ViscBndry.H>
#include <FluxBoxes.H>
#include <AMReX_MLLinOp.H>

//
// Include files for tensor solve.
//
#include <ViscBndryTensor.H>
#include <AMReX_DivVis.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MCMultiGrid.H>
#include <AMReX_MCCGSolver.H>

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

    Diffusion (amrex::Amr*               Parent,
               NavierStokesBase*  Caller,
               Diffusion*         coarser,
               int                num_state,
               amrex::FluxRegister*      Viscflux_reg,
               const amrex::Vector<int>&  _is_diffusive,
               const amrex::Vector<amrex::Real>& _visc_coef);

    ~Diffusion ();

    void echo_settings () const;

    amrex::FluxRegister* viscFluxReg ();

    amrex::Real get_scaled_abs_tol (const amrex::MultiFab& rhs,
                             amrex::Real            reduction) const;

    void diffuse_scalar (amrex::Real                   dt,
			 int                    sigma,
			 amrex::Real                   be_cn_theta,
			 const amrex::MultiFab&        rho_half,
			 int                    rho_flag,
			 amrex::MultiFab* const*       fluxn,
			 amrex::MultiFab* const*       fluxnp1,
			 int                    fluxComp,
			 amrex::MultiFab*              delta_rhs, 
                         int                    rhsComp,
			 const amrex::MultiFab*        alpha, 
                         int                    alphaComp,
			 const amrex::MultiFab* const* betan, 
			 const amrex::MultiFab* const* betanp1,
                         int                    betaComp,
			 const SolveMode&       solve_mode = ONEPASS,
                         bool                   add_old_time_divFlux = true);
    
    void diffuse_velocity (amrex::Real                   dt,
                           amrex::Real                   be_cn_theta,
                           const amrex::MultiFab&        rho_half,
                           int                    rho_flag,
                           amrex::MultiFab*              delta_rhs, 
                           int                    rhsComp,
                           const amrex::MultiFab* const* betan, 
                           const amrex::MultiFab* const* betanp1,
                           int                    betaComp);

    void diffuse_velocity (amrex::Real                   dt,
                           amrex::Real                   be_cn_theta,
                           const amrex::MultiFab&        rho_half,
                           int                    rho_flag,
                           amrex::MultiFab*              delta_rhs, 
                           const amrex::MultiFab* const* betan, 
                           const amrex::MultiFab* const* betanp1);

    void diffuse_velocity_constant_mu (amrex::Real            dt,
                                       amrex::Real            be_cn_theta,
                                       const amrex::MultiFab& rho_half,
                                       amrex::MultiFab*       delta_rhs,
                                       int             rhsComp);

    void diffuse_tensor_velocity (amrex::Real                   dt,
                                  amrex::Real                   be_cn_theta,
                                  const amrex::MultiFab&        rho_half,
                                  int                    rho_flag,
                                  amrex::MultiFab*              delta_rhs, 
                                  int                    rhsComp,
                                  const amrex::MultiFab* const* betan, 
                                  const amrex::MultiFab* const* betanp1,
                                  int                    betaComp);

    void diffuse_Vsync (amrex::MultiFab&              Vsync,
                        amrex::Real                   dt,
                        amrex::Real                   be_cn_theta,
                        const amrex::MultiFab&        rho_half,
                        int                    rho_flag,
                        const amrex::MultiFab* const* beta,
                        int                    betaComp = 0,
			bool                   update_fluxreg = true);

    void diffuse_Vsync_constant_mu (amrex::MultiFab&       Vsync,
                                    amrex::Real            dt,
                                    amrex::Real            be_cn_theta,
                                    const amrex::MultiFab& rho_half,
                                    int             rho_flag,
				    bool            update_fluxreg);

    void diffuse_tensor_Vsync (amrex::MultiFab&              Vsync,
                               amrex::Real                   dt,
                               amrex::Real                   be_cn_theta,
                               const amrex::MultiFab&        rho_half,
                               int                    rho_flag,
                               const amrex::MultiFab* const* beta,
                               int                    betaComp,
			       bool                   update_fluxreg);

    void diffuse_Ssync (amrex::MultiFab&              Ssync,
                        int                    sigma,
                        amrex::Real                   dt,
                        amrex::Real                   be_cn_theta,
                        const amrex::MultiFab&        rho_half,
                        int                    rho_flag,
			amrex::MultiFab* const*       flux,
			int                    fluxComp,
                        const amrex::MultiFab* const* beta,
                        int                    betaComp,
			const amrex::MultiFab*        alpha,
                        int                    alphaComp);

    amrex::ABecLaplacian* getViscOp (int                    src_comp,
                              amrex::Real                   a,
                              amrex::Real                   b, 
                              amrex::Real                   time,
                              ViscBndry&             visc_bndry,
                              const amrex::MultiFab&        rho_half,
                              int                    rho_flag,
                              amrex::Real*                  rhsscale,
                              const amrex::MultiFab* const* beta,
			      int                    betaComp,
                              const amrex::MultiFab*        alpha_in,
			      int                    alphaComp,
                              bool		     bndry_already_filled = false);

    amrex::ABecLaplacian* getViscOp (int                    src_comp,
                              amrex::Real                   a,
                              amrex::Real                   b,
                              const amrex::MultiFab&        rho_half,
                              int                    rho_flag,
                              amrex::Real*                  rhsscale,
                              const amrex::MultiFab* const* beta,
			      int                    betaComp,
                              const amrex::MultiFab*        alpha_in,
			      int                    alphaComp);

    void setAlpha (amrex::ABecLaplacian*  visc_op,
                   int             src_comp,
                   amrex::Real            a,
                   amrex::Real            b, 
                   amrex::Real            time,
                   const amrex::MultiFab& rho_half,
                   int                    rho_flag,
                   amrex::Real*                  rhsscale,
                   int                    dataComp,
                   const amrex::MultiFab* alpha_in);

    void setBeta (amrex::ABecLaplacian*         visc_op,
                  const amrex::MultiFab* const* beta,
                  int                    betaComp);

    amrex::DivVis* getTensorOp (amrex::Real                   a,
                         amrex::Real                   b, 
                         amrex::Real                   time, 
                         ViscBndryTensor&       visc_bndry,
                         const amrex::MultiFab&        rho_half,
                         const amrex::MultiFab* const* beta,
			 int                    betaComp);

    void getTensorBndryData (ViscBndryTensor& bndry, amrex::Real time);

    amrex::DivVis* getTensorOp (amrex::Real                   a,
                         amrex::Real                   b, 
                         const amrex::MultiFab&        rho_half,
                         const amrex::MultiFab* const* beta,
			 int                    betaComp);

    void getViscTerms (amrex::MultiFab&              visc_terms,
                       int                    src_comp,
                       int                    comp,
                       amrex::Real                   time,
                       int                    rho_flag,
                       const amrex::MultiFab* const* beta,
		       int                    betaComp);

    void getTensorViscTerms (amrex::MultiFab&              visc_terms, 
                             amrex::Real                   time,
                             const amrex::MultiFab* const* beta,
			     int                    dataComp);

    void getBndryData (ViscBndry&         bndry,
                       int                state_ind,
                       int                num_comp,
                       amrex::Real               time,
                       int                rho_flag);

    void getBndryDataGivenS (ViscBndry&         bndry,
                             amrex::MultiFab&          S,
                             amrex::MultiFab&          S_crse,
                             int                state_ind,
                             int                src_comp,
                             int                num_comp);

    void FillBoundary (amrex::BndryRegister&     bdry,
                       int                src_comp,
                       int                dest_comp,
                       int                num_comp,
                       amrex::Real               time,
                       int                rho_flag);

    void checkBetas (const amrex::MultiFab* const* beta1,
                     const amrex::MultiFab* const* beta2,
                     int&                   allthere,
                     int&                   allnull) const;

    void checkBeta (const amrex::MultiFab* const* beta,
                    int&                   allthere,
                    int&                   allnull) const;

    void checkBeta (const amrex::MultiFab* const* beta,
                    int&                   allthere) const;

    void compute_divmusi (amrex::Real                   time,
			  const amrex::MultiFab* const* beta,
			  amrex::MultiFab&              divmusi);

    void compute_divmusi (amrex::Real      time,
			  amrex::Real      mu,
			  amrex::MultiFab& divmusi);

    int maxOrder () const;
    int tensorMaxOrder () const;

    static int set_rho_flag (const DiffusionForm compDiffusionType);

    static bool are_any (const amrex::Vector<DiffusionForm>& diffusionType,
                         const DiffusionForm         testForm,
                         const int                   sComp,
                         const int                   nComp);

    static int how_many (const amrex::Vector<DiffusionForm>& diffusionType,
                         const DiffusionForm         testForm,
                         const int                   sComp,
                         const int                   nComp);

protected:

    void getTensorOp_doit (amrex::DivVis*                tensor_op,
                           amrex::Real                   a,
                           amrex::Real                   b,
                           const amrex::MultiFab&        rho,
                           const amrex::MultiFab* const* beta,
                           int                    betaComp);

    void setDomainBC (std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
                      std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc,
                      int src_comp);

    void computeAlpha (amrex::MultiFab&       alpha,
                       std::pair<amrex::Real,amrex::Real>& scalars,
                       int                    comp,
                       amrex::Real            a,
                       amrex::Real            b,
                       amrex::Real            time,
                       const amrex::MultiFab& rho,
                       int                    rho_flag, 
                       amrex::Real*           rhsscale,
                       int                    dataComp,
                       const amrex::MultiFab* alpha_in);

    void computeBeta (std::array<amrex::MultiFab,AMREX_SPACEDIM>& bcoeffs,
                      const amrex::MultiFab* const* beta, int betaComp);

    static void Finalize ();
    //
    // Data Required by Derived Classes
    //
    amrex::Amr*              parent;
    NavierStokesBase* navier_stokes;
    const amrex::BoxArray&   grids;
    const amrex::DistributionMapping& dmap;
    const int         level;
    //
    // Static data.
    //
    static bool        use_mg_precond_flag;
    static int         use_mlmg_solver;
    static int         use_cg_solve;
    static int         scale_abec;
    static amrex::Vector<int>  is_diffusive;    // Does variable diffuse?
    static amrex::Vector<amrex::Real> visc_coef;       // Const coef viscosity terms
    static int         verbose;
    static amrex::Real        visc_tol;

private:
    //
    // The data.
    //
    Diffusion*    coarser;
    Diffusion*    finer;
    int           NUM_STATE;
    amrex::IntVect       crse_ratio;
    amrex::FluxRegister* viscflux_reg;
    //
    // Static data.
    //
    static int         do_reflux;
    static int         use_tensor_cg_solve;
    static int         max_order;
    static int         tensor_max_order;
};

#endif

