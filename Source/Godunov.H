
#ifndef _GODUNOV_CLASS_H_
#define _GODUNOV_CLASS_H_ 

//
// A useful enumeration of the forms for the advective terms
//
enum AdvectionForm {Conservative = 0, NonConservative};
enum AdvectionScheme {PRE_MAC, FPU, BDS};

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
    Godunov( int max_size = 512 );

    // destruct the advection object
    ~Godunov();

    // set up the Godunov amrex::Box
    void BuildWorkSpace (const amrex::Box& grd, const amrex::Real* dx, amrex::Real dt);

    void AllocEdgeBoxes (const amrex::Box& grd, D_DECL(amrex::FArrayBox& xflux,
                                                        amrex::FArrayBox& yflux,
                                                        amrex::FArrayBox& zflux));

    void ComputeTransverVelocities (const amrex::Box& grd, const amrex::Real* dx, amrex::Real dt,
                                    D_DECL(const int* ubc,const int* vbc,const int* wbc),
                                    D_DECL(const amrex::FArrayBox& U,const amrex::FArrayBox& V,const amrex::FArrayBox& W),
                                    D_DECL(int Ucomp,int Vcomp,int Wcomp),
                                    const amrex::FArrayBox& tforces, int Tcomp);

    void Setup (const amrex::Box& grd, const amrex::Real* dx, amrex::Real dt, int velpred, 
                amrex::FArrayBox& xflux, const int* ubc,
                amrex::FArrayBox& yflux, const int* vbc,
#if (BL_SPACEDIM == 3 )
                amrex::FArrayBox& zflux, const int* wbc,
#endif
                const amrex::FArrayBox& U, const amrex::FArrayBox& Rho, 
                const amrex::FArrayBox& tforces);

    void Setup (const amrex::Box& grd, const amrex::Real* dx, amrex::Real dt, int velpred,
                D_DECL(amrex::FArrayBox& xflux,amrex::FArrayBox& yflux,amrex::FArrayBox& zflux),
                D_DECL(const int* ubc,const int* vbc,const int* wbc),
                D_DECL(const amrex::FArrayBox& U,const amrex::FArrayBox& V,const amrex::FArrayBox& W),
                D_DECL(int Ucomp,int Vcomp,int Wcomp),
                const amrex::FArrayBox& tforces, int Tcomp);

    // =============================================================
    // Data access functions
    // ==============================================================

    int useForcesInTrans() {return use_forces_in_trans;}
    
    // =============================================================
    // Advection functions follow
    // =============================================================

    void edge_states (const amrex::Box& grd, const amrex::Real* dx, amrex::Real dt, int velpred,
                      const amrex::FArrayBox& uedge, amrex::FArrayBox& stx,
                      const amrex::FArrayBox& vedge, amrex::FArrayBox& sty,
#if (BL_SPACEDIM == 3 )
                      const amrex::FArrayBox& wedge, amrex::FArrayBox& stz,
#endif
                      const amrex::FArrayBox& U, const amrex::FArrayBox& S, 
                      const amrex::FArrayBox& tforces, const amrex::FArrayBox& divu,
                      int fab_ind, int state_ind, const int *bc,
                      int use_conserv, AdvectionScheme which_scheme);

    void edge_states (const amrex::Box& grd, const amrex::Real* dx, amrex::Real dt, int velpred,
                      D_DECL(const amrex::FArrayBox &uedge, const amrex::FArrayBox &vedge, const amrex::FArrayBox &wedge),
                      D_DECL(int mCompX,       int mCompY,       int mCompZ),
                      D_DECL(amrex::FArrayBox &stx,   amrex::FArrayBox &sty,   amrex::FArrayBox &stz  ),
                      D_DECL(int eCompX,       int eCompY,       int eCompZ),
                      D_DECL(const amrex::FArrayBox& U, const amrex::FArrayBox &V, const amrex::FArrayBox &W),
                      D_DECL(int Ucomp,        int Vcomp,        int Wcomp),
                      const amrex::FArrayBox &S, int Scomp, const amrex::FArrayBox &tforces, int Tcomp, const amrex::FArrayBox& divu, int Dcomp,
                      int state_ind, const int *bc, int use_conserv, AdvectionScheme which_scheme);

    void edge_states_orig( const amrex::Box &grd, const amrex::Real *dx, amrex::Real dt, int velpred,
                           D_DECL(const amrex::FArrayBox &uedge, const amrex::FArrayBox &vedge, const amrex::FArrayBox &wedge),
                           D_DECL(int mCompX,       int mCompY,       int mCompZ),
                           D_DECL(amrex::FArrayBox &stx,   amrex::FArrayBox &sty,   amrex::FArrayBox &stz  ),
                           D_DECL(int eCompX,       int eCompY,       int eCompZ),
                           D_DECL(const amrex::FArrayBox& U, const amrex::FArrayBox &V, const amrex::FArrayBox &W),
                           D_DECL(int Ucomp,        int Vcomp,        int Wcomp),
                           const amrex::FArrayBox &S, int Scomp, const amrex::FArrayBox &tforces, int Tcomp,
                           int state_ind, const int *bc);

    void edge_states_fpu( const amrex::Box &grd, const amrex::Real *dx, amrex::Real dt,
                          D_DECL(const amrex::FArrayBox &uedge, const amrex::FArrayBox &vedge, const amrex::FArrayBox &wedge),
                          D_DECL(int mCompX,       int mCompY,       int mCompZ),
                          D_DECL(amrex::FArrayBox &stx,   amrex::FArrayBox &sty,   amrex::FArrayBox &stz  ),
                          D_DECL(int eCompX,       int eCompY,       int eCompZ),
                          const amrex::FArrayBox &S, int Scomp, const amrex::FArrayBox &tforces, int Tcomp, const amrex::FArrayBox& divu, int Dcomp,
                          int state_ind, const int *bc, int use_conserv);
    
    void edge_states_bds( const amrex::Box &grd, const amrex::Real *dx, amrex::Real dt,
                          D_DECL(const amrex::FArrayBox &uedge, const amrex::FArrayBox &vedge, const amrex::FArrayBox &wedge),
                          D_DECL(int mCompX,       int mCompY,       int mCompZ),
                          D_DECL(amrex::FArrayBox &stx,   amrex::FArrayBox &sty,   amrex::FArrayBox &stz  ),
                          D_DECL(int eCompX,       int eCompY,       int eCompZ),
                          const amrex::FArrayBox &S, int Scomp, const amrex::FArrayBox &tforces, int Tcomp, const amrex::FArrayBox& divu, int Dcomp,
                          int state_ind, const int *bc, int use_conserv);
    
    // compute divergent mac velocities
    void ComputeUmac( const amrex::Box &grd, const amrex::Real *dx, amrex::Real dt, 
                      amrex::FArrayBox &umac, const int *ubc, 
                      amrex::FArrayBox &vmac, const int *vbc, 
#if (BL_SPACEDIM == 3 )
                      amrex::FArrayBox &wmac, const int *wbc, 
#endif
                      amrex::FArrayBox &U, amrex::FArrayBox &tforces );
    
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

    // test the cell centered Courant number 
    amrex::Real test_u_rho( amrex::FArrayBox &U, amrex::FArrayBox &Rho, const amrex::Box &grd,
                     const amrex::Real *dx, const amrex::Real dt,
                     const amrex::Real *u_max );

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

    // correct the 1st order RK to 2nd order via
    // psi^n+1 = psi^* + (dt/2)*(tf^* - tf^n)
    void Correct_tf( const amrex::FArrayBox &Sstar, amrex::FArrayBox &Snew,
                     int start_ind, int num_comp, 
                     const amrex::FArrayBox &tfstar, const amrex::FArrayBox &tfn,
                     int tf_ind,
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
    
    // manipulate 1d scratch space
    void SetBogusScratch();
    void SetScratch( int max_size );

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

