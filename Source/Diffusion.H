
#ifndef IAMR_Diffusion_H_
#define IAMR_Diffusion_H_


#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_LevelBld.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_ErrorList.H>
#include <AMReX_FluxRegister.H>
#include <FluxBoxes.H>

#include <AMReX_MLMG.H>
#ifdef AMREX_USE_EB
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLEBTensorOp.H>
#else
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLTensorOp.H>
#endif


//
// Useful enumeration of the different forms of the diffusion terms
//
enum DiffusionForm { RhoInverse_Laplacian_S, Laplacian_SoverRho, Laplacian_S };

class NavierStokesBase;

class Diffusion
{
public:

    Diffusion ();

    Diffusion (amrex::Amr*               Parent,
               NavierStokesBase*         Caller,
               Diffusion*                coarser,
               int                       num_state,
               amrex::FluxRegister*      Viscflux_reg,
               const amrex::Vector<int>&  _is_diffusive);

    void echo_settings () const;

    amrex::FluxRegister* viscFluxReg ();

    static amrex::Real get_scaled_abs_tol (const amrex::MultiFab& rhs,
                                           amrex::Real            reduction);

    void diffuse_scalar (const amrex::Vector<amrex::MultiFab*>& S_old,
                         const amrex::Vector<amrex::MultiFab*>& Rho_old,
                         amrex::Vector<amrex::MultiFab*>&       S_new,
                         const amrex::Vector<amrex::MultiFab*>& Rho_new,
                         int                                    S_comp,
                         int                                    nComp,
                         int                                    Rho_comp,
                         amrex::Real                            prev_time,
                         amrex::Real                            curr_time,
                         amrex::Real                            be_cn_theta,
                         const amrex::MultiFab&                 rho_half,
                         int                                    rho_flag,
                         amrex::MultiFab* const*                fluxn,
                         amrex::MultiFab* const*                fluxnp1,
                         int                                    fluxComp,
                         amrex::MultiFab*                       delta_rhs,
                         int                                    rhsComp,
                         const amrex::MultiFab*                 alpha,
                         int                                    alpha_in_comp,
                         const amrex::MultiFab* const*          betan,
                         const amrex::MultiFab* const*          betanp1,
                         int                                    betaComp,
                         const amrex::IntVect&                  cratio,
                         const amrex::BCRec&                    bc,
                         const amrex::Geometry&                 geom,
                         bool                                   add_old_time_divFlux = true,
                         const amrex::Vector<int>&              is_diffusive = amrex::Vector<int>());

   void diffuse_velocity (amrex::Real                   dt,
                          amrex::Real                   be_cn_theta,
                          const amrex::MultiFab&        rho_half,
                          int                           rho_flag,
                          amrex::MultiFab*              delta_rhs,
                          const amrex::MultiFab* const* betan,
                          const amrex::MultiFab*        betanCC,
                          const amrex::MultiFab* const* betanp1,
                          const amrex::MultiFab*        betanp1CC);

    void diffuse_velocity (amrex::Real                   dt,
                           amrex::Real                   be_cn_theta,
                           const amrex::MultiFab&        rho_half,
                           int                           rho_flag,
                           amrex::MultiFab*              delta_rhs,
                           int                           rhsComp,
                           const amrex::MultiFab* const* betan,
                           const amrex::MultiFab*        betanCC,
                           const amrex::MultiFab* const* betanp1,
                           const amrex::MultiFab*        betanp1CC,
                           int                           betaComp);

    void diffuse_tensor_velocity (amrex::Real            dt,
                                  amrex::Real            be_cn_theta,
                                  const amrex::MultiFab& rho_half,
                                  int                    rho_flag,
                                  amrex::MultiFab*       delta_rhs,
                                  int                    rhsComp,
                                  const amrex::MultiFab* const* betan,
                                  const amrex::MultiFab*        betanCC,
                                  const amrex::MultiFab* const* betanp1,
                                  const amrex::MultiFab*        betanp1CC,
                                  int                    betaComp);

    void diffuse_Vsync (amrex::MultiFab&              Vsync,
                        amrex::Real                   dt,
                        amrex::Real                   be_cn_theta,
                        const amrex::MultiFab&        rho_half,
                        int                           rho_flag,
                        const amrex::MultiFab*        const* beta,
                        int                           betaComp = 0,
                        bool                          update_fluxreg = true);

    void diffuse_tensor_Vsync (amrex::MultiFab&              Vsync,
                               amrex::Real                   dt,
                               amrex::Real                   be_cn_theta,
                               const amrex::MultiFab&        rho_half,
                               int                           rho_flag,
                               const amrex::MultiFab* const* beta,
                               int                           betaComp,
                               bool                          update_fluxreg);


    void diffuse_Ssync (amrex::MultiFab&              Ssync,
                        int                           sigma,
                        int                           nComp,
                        amrex::Real                   dt,
                        amrex::Real                   be_cn_theta,
                        const amrex::MultiFab&        rho_half,
                        int                           rho_flag,
                        amrex::MultiFab* const*       flux,
                        int                           fluxComp,
                        const amrex::MultiFab* const* beta,
                        int                           betaComp,
                        const amrex::MultiFab*        alpha,
                        int                           alphaComp);


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
                             const amrex::MultiFab*        betaCC,
                             int                           betaComp);

    void FillBoundary (amrex::BndryRegister& bdry,
                       int                   state_ind,
                       int                   dest_comp,
                       int                   num_comp,
                       amrex::Real           time,
                       int                   rho_flag);

    static void checkBeta (const amrex::MultiFab* const* beta,
                           int&                          allthere,
                           int&                          allnull);

    static void checkBeta (const amrex::MultiFab* const* beta,
                           int&                          allthere);

    [[nodiscard]] static int maxOrder ();
    [[nodiscard]] static int tensorMaxOrder ();

    static int set_rho_flag (DiffusionForm compDiffusionType);

    static void setDomainBC (std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
                             std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc,
                             const amrex::BCRec&                            bc);

    static void computeAlpha (amrex::MultiFab&                    alpha,
                              std::pair<amrex::Real,amrex::Real>& scalars,
                              amrex::Real                         a,
                              amrex::Real                         b,
                              amrex::Real*                        rhsscale,
                              const amrex::MultiFab*              alpha_in,
                              int                                 alpha_in_comp,
                              int                                 rho_flag,
                              const amrex::MultiFab*              rho,
                              int                                 rho_comp);

#ifdef AMREX_USE_EB
    static void setBeta(amrex::MLEBABecLap&           op,
                        const amrex::MultiFab* const* beta,
                        int                           betaComp,
                        int                           nComp = 1);

    static void setViscosity(amrex::MLEBTensorOp&          tensorop,
                      const amrex::MultiFab* const* beta,
                      int                           betaComp,
                      const amrex::MultiFab&        beta_cc);
#else
    static void setBeta(amrex::MLABecLaplacian&       op,
                        const amrex::MultiFab* const* beta,
                        int                           betaComp,
                        int                           nComp = 1);

    static void setViscosity(amrex::MLTensorOp&            tensorop,
                             const amrex::MultiFab* const* beta,
                             int                           betaComp);
#endif

    static void computeExtensiveFluxes(amrex::MLMG&            a_mg,
                                       amrex::MultiFab&        Soln,
                                       amrex::MultiFab* const* flux,
                                       int                     fluxComp,
                                       int                     ncomp,
                                       amrex::MultiFab*        area,
                                       amrex::Real             fac );

protected:

    void setDomainBC (std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
                      std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc,
                      int src_comp);


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
    static int         scale_abec;
    static amrex::Vector<int>  is_diffusive;    // Does variable diffuse?
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
    static int         max_order;
    static int         tensor_max_order;
};

#endif
