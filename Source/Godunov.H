
#ifndef _GODUNOV_CLASS_H_
#define _GODUNOV_CLASS_H_

#ifdef AMREX_USE_EB
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_BCRec.H>
#endif

//
// A useful enumeration of the forms for the advective terms
//
enum AdvectionForm {Conservative = 0, NonConservative};
enum AdvectionScheme {PRE_MAC, FPU, BDS};

#define COVERED_VAL 1.0e40

//===========================================================
// This header file defines the Godunov driver class.  Its
// interface simplifies the calling sequence for the Godunov
// box
//===========================================================

class Godunov {
public:      // public access functions

    // =============================================================
    // Setup functions follow
    // =============================================================

    // construct the advection object
    Godunov(int max_size = 512);


    // destruct the advection object
    ~Godunov();

    void ExtrapVelToFaces (const amrex::Box&  box,
                           const amrex_real*  dx,
                           amrex_real         dt,
                           D_DECL(      amrex::FArrayBox&  umac,       amrex::FArrayBox&  vmac,       amrex::FArrayBox&  wmac),
                           D_DECL(const amrex::Vector<int>& ubc, const amrex::Vector<int>& vbc, const amrex::Vector<int>& wbc),
                           const amrex::FArrayBox&  U,
                           amrex::FArrayBox&  tforces);

    void AdvectScalars(const amrex::Box&  box,
                       const amrex::Real* dx,
                       amrex::Real        dt,
                       D_DECL(const amrex::FArrayBox&   Ax, const amrex::FArrayBox&   Ay, const amrex::FArrayBox&   Az),
                       D_DECL(const amrex::FArrayBox& umac, const amrex::FArrayBox& vmac, const amrex::FArrayBox& wmac),   int vcomp,
                       D_DECL(      amrex::FArrayBox& xflx,       amrex::FArrayBox& yflx,       amrex::FArrayBox& zflx),   int flxcomp,
                       D_DECL(      amrex::FArrayBox& xstate,     amrex::FArrayBox& ystate,     amrex::FArrayBox& zstate), int ecomp,
                       const amrex::FArrayBox& Sfab,   int first_scalar, int num_scalars,
                       const amrex::FArrayBox& Forces, int fcomp,
                       const amrex::FArrayBox& Divu,   int ducomp,
                       amrex::FArrayBox& aofs,         int state_ind,
                       const amrex::Vector<AdvectionForm>& advectionType, const amrex::Vector<int>& state_bc,
                       AdvectionScheme adv_scheme, const amrex::FArrayBox& V);

    // =============================================================
    // Data access functions
    // ==============================================================

    int useForcesInTrans() {return use_forces_in_trans;}

    // =============================================================
    // Advection functions follow
    // =============================================================

    // advect a state component
    void AdvectState( const amrex::Box &grd, const amrex::Real *dx, amrex::Real dt,
                      amrex::FArrayBox &areax, amrex::FArrayBox &uedge, amrex::FArrayBox &xflux,
                      amrex::FArrayBox &areay, amrex::FArrayBox &vedge, amrex::FArrayBox &yflux,
#if (BL_SPACEDIM == 3 )
                      amrex::FArrayBox &areaz, amrex::FArrayBox &wedge, amrex::FArrayBox &zflux,
#endif
                      amrex::FArrayBox &U,
                      amrex::FArrayBox &S,
                      amrex::FArrayBox &tforces, amrex::FArrayBox& divu,
                      int fab_ind,
                      amrex::FArrayBox &aofs,
                      int aofs_ind,
                      int iconserv, int state_ind, const int *bc,
                      AdvectionScheme scheme, amrex::FArrayBox &vol );

    // compute the advective derivative from fluxes
    void ComputeAofs (const amrex::Box& grd,
                      const amrex::FArrayBox& areax, const amrex::FArrayBox& uedge, const amrex::FArrayBox& xflux,
                      const amrex::FArrayBox& areay, const amrex::FArrayBox& vedge, const amrex::FArrayBox& yflux,
#if (BL_SPACEDIM == 3 )
                      const amrex::FArrayBox& areaz, const amrex::FArrayBox& wedge, const amrex::FArrayBox& zflux,
#endif
                      const amrex::FArrayBox& vol,
                      amrex::FArrayBox& aofs, int aofs_ind, int iconserv) const;

    void ComputeAofs (const amrex::Box& grd,
                      D_DECL(const amrex::FArrayBox& areax,
                             const amrex::FArrayBox& areay,
                             const amrex::FArrayBox& areaz),
                      D_DECL(int axcomp,int aycomp,int azcomp),
                      D_DECL(const amrex::FArrayBox& uedge,
                             const amrex::FArrayBox& vedge,
                             const amrex::FArrayBox& wedge),
                      D_DECL(int ucomp, int vcomp, int wcomp),
                      D_DECL(const amrex::FArrayBox& xflux,
                             const amrex::FArrayBox& yflux,
                             const amrex::FArrayBox& zflux),
                      D_DECL(int fxcomp,int fycomp,int fzcomp),
                      const amrex::FArrayBox& vol, int volcomp,
                      amrex::FArrayBox& aofs,int acomp, int iconserv ) const;

    // sync advect a state component
    void SyncAdvect( const amrex::Box &grd, const amrex::Real *dx, amrex::Real dt, int level,
                     const amrex::FArrayBox &areax, amrex::FArrayBox &uedge,
                     amrex::FArrayBox &ucorr, amrex::FArrayBox &xflux,
                     const amrex::FArrayBox &areay, amrex::FArrayBox &vedge,
                     amrex::FArrayBox &vcorr, amrex::FArrayBox &yflux,
#if (BL_SPACEDIM == 3 )
                     const amrex::FArrayBox &areaz, amrex::FArrayBox &wedge,
                     amrex::FArrayBox &wcorr, amrex::FArrayBox &zflux,
#endif
                     amrex::FArrayBox &U, amrex::FArrayBox &S,
                     amrex::FArrayBox &tforces, amrex::FArrayBox& divu, int ind,
                     amrex::FArrayBox &sync,                     int sync_ind,
                     int iconserv, int state_ind, const int *bc,
                     AdvectionScheme scheme, const amrex::FArrayBox &vol );

    // compute the advective derivative of corrective fluxes for the mac sync
    void ComputeSyncAofs( const amrex::Box &grd,
                          const amrex::FArrayBox &areax, amrex::FArrayBox &ucorr, amrex::FArrayBox &xflux,
                          const amrex::FArrayBox &areay, amrex::FArrayBox &vcorr, amrex::FArrayBox &yflux,
#if (BL_SPACEDIM == 3 )
                          const amrex::FArrayBox &areaz, amrex::FArrayBox &wcorr, amrex::FArrayBox &zflux,
#endif
                          const amrex::FArrayBox &vol,
                          amrex::FArrayBox &sync,
                          int sync_ind, int iconserv );

    // correct a conservatively-advected scalar for under-over shoots
    void ConservativeScalMinMax( amrex::FArrayBox &Sold, amrex::FArrayBox &Snew,
                                 int ind_old_s, int ind_old_rho,
                                 int ind_new_s, int ind_new_rho,
                                 const int *bc, const amrex::Box &grd );

    // correct a convectively-advected scalar for under-over shoots
    void ConvectiveScalMinMax( amrex::FArrayBox &Sold, amrex::FArrayBox &Snew,
                               int ind_old, int ind_new,
                               const int *bc, const amrex::Box &grd );

    // =============================================================
    // Diagnostic functions follow
    // =============================================================

    // estimate the timestep from the cell-centered quantities
    amrex::Real estdt( amrex::FArrayBox &U, amrex::FArrayBox &tforces, amrex::FArrayBox &Rho,
                const amrex::Box &grd, const amrex::Real *dx, amrex::Real cfl, amrex::Real *u_max );

	// find the largest change in velocity magnitude since last iteration
    amrex::Real maxchng_velmag( amrex::FArrayBox &U_old, amrex::FArrayBox &U_new, const amrex::Box &grd );

    // test the mac edge velocity Courant numbers
    amrex::Real test_umac_rho( amrex::FArrayBox &umac,
                        amrex::FArrayBox &vmac,
#if (BL_SPACEDIM == 3 )
                        amrex::FArrayBox &wmac,
#endif
                        amrex::FArrayBox &Rho,
                        const amrex::Box &grd, const amrex::Real *dx, const amrex::Real dt,
                        const amrex::Real *u_max );

    static bool are_any(const amrex::Vector<AdvectionForm>& advectionType,
                        const AdvectionForm         testForm,
                        const int                   sComp,
                        const int                   nComp);

    static int how_many(const amrex::Vector<AdvectionForm>& advectionType,
                        const AdvectionForm         testForm,
                        const int                   sComp,
                        const int                   nComp);

    static int hypgrow ();

    // =============================================================
    // Source term functions follow
    // =============================================================

    // compute the update rule, this is useful for 1st order RK
    // psi^n+1 = psi^n + dt*(tf^n)
    void Add_tf( const amrex::FArrayBox &Sold,
                 amrex::FArrayBox &Snew,    int start_ind, int num_comp,
                 const amrex::FArrayBox &tforces, int tf_ind,
                 const amrex::Box &grd,     amrex::Real dt ) const;

    // compute the update rule
    // psi^n+1 = psi^n - dt*aofs + dt*tforces
    void Add_aofs_tf( const amrex::FArrayBox &Sold,
                      amrex::FArrayBox &Snew,    int start_ind, int num_comp,
                      const amrex::FArrayBox &Aofs,    int aofs_ind,
                      const amrex::FArrayBox &tforces, int tf_ind,
                      const amrex::Box &grd,     amrex::Real dt ) const;

    // compute the update rule for velocities or momenta
    // psi^n+1 = psi^n - dt*aofs - dt*gp/rho + dt*tforces
    void Add_aofs_tf_gp( const amrex::FArrayBox &Uold, amrex::FArrayBox &Unew,
                         const amrex::FArrayBox &Aofs, const amrex::FArrayBox &tforces,
                         const amrex::FArrayBox &gp,   const amrex::FArrayBox &rho,
                         const amrex::Box &grd,  amrex::Real dt) const;

    // compute total source term for velocities
    // tforces = (tforces + visc - gp)/rho
    void Sum_tf_gp_visc (amrex::FArrayBox&       tforces, int Tcomp,
                         const amrex::FArrayBox& visc,    int Vcomp,
                         const amrex::FArrayBox& gp,      int Gcomp,
                         const amrex::FArrayBox& rho,     int Rcomp) const;

    void Sum_tf_gp_visc( amrex::FArrayBox&       tforces,
                         const amrex::FArrayBox& visc,
                         const amrex::FArrayBox& gp,
                         const amrex::FArrayBox& Rho) const;

    // iconserv==1   => tforces = tforces + visc - divU*S
    // iconserv==0   => tforces = (tforces+visc)/rho
    void Sum_tf_divu_visc (const amrex::FArrayBox& S,
                           amrex::FArrayBox& tforces,
                           int s_ind, int num_comp,
                           const amrex::FArrayBox& visc, int v_ind,
                           const amrex::FArrayBox& divu,
                           const amrex::FArrayBox& rho,
                           int iconserv) const;

    void Sum_tf_divu_visc (const amrex::FArrayBox& S, int s_ind,
                           amrex::FArrayBox& tforces, int t_ind,
                           int        num_comp,
                           const amrex::FArrayBox& visc,    int v_ind,
                           const amrex::FArrayBox& divu,    int d_ind,
                           const amrex::FArrayBox& rho,     int tr_ind,
                           int        iconserv) const;

    // compute total source term for velocities
    // tforces = (tforces - gp)/rho
    void Sum_tf_gp (amrex::FArrayBox& tforces, int Tcomp,
                    const amrex::FArrayBox& gp, int Gcomp,
                    const amrex::FArrayBox& rho, int Rcomp) const;

    // iconserv==1   => tforces = tforces - divU*S
    // iconserv==0   => tforces = (tforces)/rho
    void Sum_tf_divu (const amrex::FArrayBox& S,       int s_ind,
                      amrex::FArrayBox& tforces, int t_ind,
                      int        num_comp,
                      const amrex::FArrayBox& divu,    int d_ind,
                      const amrex::FArrayBox& rho,     int r_ind,
                      int        iconserv) const;

  //
  // EB functions
  //
#ifdef AMREX_USE_EB

    void ComputeSlopes( amrex::MultiFab&  a_Sborder,
                        D_DECL(amrex::MultiFab& a_xslopes,
                               amrex::MultiFab& a_yslopes,
                               amrex::MultiFab& a_zslopes),
                        const  amrex::Vector<amrex::BCRec>& a_bcs,
                        int a_comp,
                        int a_ncomp,
                        const amrex::Box& a_domain);


  // void ComputeVelocitySlopes(const amrex::MFIter&  mfi,
  //                            const amrex::MultiFab& Sborder,
  //                            D_DECL(amrex::MultiFab& m_xslopes,
  //                                   amrex::MultiFab& m_yslopes,
  //                                   amrex::MultiFab& m_zslopes),
  //                            D_DECL(const amrex::Vector<int>& a_ubc,
  //                                   const amrex::Vector<int>& a_vbc,
  //                                   const amrex::Vector<int>& a_wbc),
  //                            const amrex::Box& domain);

  // void ComputeScalarSlopes(const amrex::MFIter& mfi,
  //                          const amrex::MultiFab& Sborder, const int& ncomp,
  //                          D_DECL(amrex::MultiFab& xslopes,
  //                                 amrex::MultiFab& yslopes,
  //                                 amrex::MultiFab& zslopes),
  //                          D_DECL(const amrex::Vector<int>& ubc,
  //                                 const amrex::Vector<int>& vbc,
  //                                 const amrex::Vector<int>& wbc),
  //                                 const amrex::Box& domain);

  // void ExtrapVelToFaces (const amrex::MFIter&  mfi,
  //                        D_DECL(amrex::FArrayBox&  umac, amrex::FArrayBox&  vmac,
  //                               amrex::FArrayBox&  wmac),
  //                        D_DECL(const amrex::Vector<int>& ubc,
  //                               const amrex::Vector<int>& vbc,
  //                               const amrex::Vector<int>& wbc),
  //                        D_DECL(amrex::MultiFab& m_xslopes,
  //                               amrex::MultiFab& m_yslopes,
  //                               amrex::MultiFab& m_zslopes),
  //                        const amrex::FArrayBox&  U,
  //                        const amrex::FArrayBox&  tforces,
  //                        D_DECL(const amrex::MultiCutFab& x_areafrac,
  //       			                        const amrex::MultiCutFab& y_areafrac,
  //       			                        const amrex::MultiCutFab& z_areafrac),
  //                        D_DECL(const amrex::MultiCutFab& x_facecent,
  //       			                        const amrex::MultiCutFab& y_facecent,
  //       			                        const amrex::MultiCutFab& z_facecent),
  //                        const amrex::Box& domain);

    void ExtrapVelToFaces ( const amrex::MultiFab&  a_vel,
                            D_DECL(amrex::MultiFab& a_umac,
                                   amrex::MultiFab& a_vmac,
                                   amrex::MultiFab& a_wmac ),
                            D_DECL(amrex::MultiFab& a_xslopes,
                                   amrex::MultiFab& a_yslopes,
                                   amrex::MultiFab& a_zslopes),
                            const amrex::Geometry&  a_geom,
                            const amrex::Vector<amrex::BCRec>& a_bcs );

    void ComputeFluxes(  D_DECL(amrex::MultiFab& a_fx,
                                amrex::MultiFab& a_fy,
                                amrex::MultiFab& a_fz),
                         amrex::MultiFab& a_state,
                         const int a_comp,
                         const int a_ncomp,
                         D_DECL( const amrex::MultiFab& a_xsl,
                                 const amrex::MultiFab& a_ysl,
                                 const amrex::MultiFab& a_zsl),
                         const int a_sl_comp,
                         D_DECL( const amrex::MultiFab& u_mac,
                                 const amrex::MultiFab& v_mac,
                                 const amrex::MultiFab& w_mac),
                         const amrex::Geometry& a_geom,
                         const amrex::Vector<amrex::BCRec>& a_bcs);


    void ComputeConvectiveTerm (amrex::MultiFab& a_state,
                                const int a_state_comp,
                                amrex::MultiFab& a_conv,
                                const int a_conv_comp,
                                const int a_ncomp,
                                D_DECL(amrex::MultiFab& a_fx,
                                       amrex::MultiFab& a_fy,
                                       amrex::MultiFab& a_fz),
                                D_DECL(const amrex::MultiFab& a_umac,
                                       const amrex::MultiFab& a_vmac,
                                       const amrex::MultiFab& a_wmac),
                                D_DECL(const amrex::MultiFab& a_xsl,
                                       const amrex::MultiFab& a_ysl,
                                       const amrex::MultiFab& a_zsl),
                                const int a_sl_comp,
                                const amrex::Vector<amrex::BCRec>& a_bcs,
                                const amrex::Geometry& a_geom );


    void AdvectVel (const amrex::MFIter&  mfi,
                    const amrex::MultiFab& U, amrex::MultiFab& aofs,
                    D_DECL(const amrex::FArrayBox& uedge,
                           const amrex::FArrayBox& vedge,
                           const amrex::FArrayBox& wedge),
                    D_DECL(      amrex::FArrayBox& xflx,
                                 amrex::FArrayBox& yflx,
                                 amrex::FArrayBox& zflx),
                    D_DECL(      amrex::FArrayBox& xstate,
                                 amrex::FArrayBox& ystate,
                                 amrex::FArrayBox& zstate),
                    D_DECL(const amrex::Vector<int>& ubc,
                           const amrex::Vector<int>& vbc,
                           const amrex::Vector<int>& wbc),
                    D_DECL(amrex::MultiFab& m_xslopes,
                           amrex::MultiFab& m_yslopes,
                           amrex::MultiFab& m_zslopes),
                    const amrex::MultiFab& volfrac,
                    const amrex::MultiCutFab& bndrycent,
                    D_DECL(const amrex::MultiCutFab& x_areafrac,
                           const amrex::MultiCutFab& y_areafrac,
                           const amrex::MultiCutFab& z_areafrac),
                    D_DECL(const amrex::MultiCutFab& x_facecent,
                           const amrex::MultiCutFab& y_facecent,
                           const amrex::MultiCutFab& z_facecent),
                    const amrex::Box& domain,
                    const amrex::Real* dx,
                    int nghost);

  void AdvectScalars_EB (const amrex::MFIter&  mfi,
                         const amrex::MultiFab& S, int first_scal, int scomp,
                         amrex::MultiFab& aofs, int acomp, int flxcomp,
                         D_DECL(const amrex::MultiFab& xslopes,
                         const amrex::MultiFab& yslopes,
                         const amrex::MultiFab& zslopes),
                         D_DECL(const amrex::FArrayBox& uedge,
                                const amrex::FArrayBox& vedge,
                                const amrex::FArrayBox& wedge),
                         D_DECL(      amrex::FArrayBox& xflx,
                                      amrex::FArrayBox& yflx,
                                      amrex::FArrayBox& zflx),
                         D_DECL(      amrex::FArrayBox& xstate,
                                      amrex::FArrayBox& ystate,
                                      amrex::FArrayBox& zstate),
                         const amrex::MultiFab& volfrac,
                         const amrex::MultiCutFab& bndrycent,
                         D_DECL(const amrex::MultiCutFab& x_areafrac,
				                        const amrex::MultiCutFab& y_areafrac,
				                        const amrex::MultiCutFab& z_areafrac),
                         D_DECL(const amrex::MultiCutFab& x_facecent,
				                        const amrex::MultiCutFab& y_facecent,
				                        const amrex::MultiCutFab& z_facecent),
                                const amrex::Vector<int>& state_bc,
		                     const amrex::Box& domain,
		                     const amrex::Real* dx,
		                     int nghost, int known_edgestate);

#endif

  //================ for debugging purposes this is made public
    // internal FABS
    amrex::FArrayBox work; // work space for the Godunov Box
    amrex::FArrayBox uad;  // work space for advective transverse velocities
    amrex::FArrayBox vad;
#if (BL_SPACEDIM == 3 )
    amrex::FArrayBox wad;
#endif

    // More workspace for modified advection scheme allowing for ppm option
    //  (pulled in to avoid f90 allocates)
    amrex::FArrayBox smp, dsvl, I;
    amrex::FArrayBox D_DECL(sedgex, sedgey, sedgez);

    // 1D arrays used in computing slopes and edges states
    amrex::Vector<amrex::Real> stxlo;
    amrex::Vector<amrex::Real> stxhi;
    amrex::Vector<amrex::Real> slxscr;
    amrex::Vector<amrex::Real> stylo;
    amrex::Vector<amrex::Real> styhi;
    amrex::Vector<amrex::Real> slyscr;
#if (BL_SPACEDIM == 3 )
    amrex::Vector<amrex::Real> stzlo;
    amrex::Vector<amrex::Real> stzhi;
    amrex::Vector<amrex::Real> slzscr;
#endif

protected:   // protected work functions

    static void Initialize ();
    static void Finalize ();

protected:   // protected data members

    // sizing boxes
    amrex::Box work_bx;
    amrex::Box xflux_bx;
    amrex::Box yflux_bx;
#if (BL_SPACEDIM == 3 )
    amrex::Box zflux_bx;
#endif

// flags

    static int verbose;
    static int use_forces_in_trans;

// algorithm control parameters

    static int  slope_order;
    static int  ppm_type;
    static int  corner_couple;

};

#endif
