
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
    Godunov();
    
    // size 1d work space explicitly
    Godunov( int max_size );

    // destruct the advection object
    ~Godunov();

    // set up the Godunov Box
    void BuildWorkSpace (const Box& grd, const Real* dx, Real dt);

    void AllocEdgeBoxes (const Box& grd, D_DECL(FArrayBox& xflux,
                                                        FArrayBox& yflux,
                                                        FArrayBox& zflux));

    void ComputeTransverVelocities (const Box& grd, const Real* dx, Real dt,
                                    D_DECL(const int* ubc,const int* vbc,const int* wbc),
                                    D_DECL(const FArrayBox& U,const FArrayBox& V,const FArrayBox& W),
                                    D_DECL(int Ucomp,int Vcomp,int Wcomp),
                                    const FArrayBox& tforces, int Tcomp);

    void Setup (const Box& grd, const Real* dx, Real dt, int velpred, 
                FArrayBox& xflux, const int* ubc,
                FArrayBox& yflux, const int* vbc,
#if (BL_SPACEDIM == 3 )
                FArrayBox& zflux, const int* wbc,
#endif
                const FArrayBox& U, const FArrayBox& Rho, 
                const FArrayBox& tforces);

    void Setup (const Box& grd, const Real* dx, Real dt, int velpred,
                D_DECL(FArrayBox& xflux,FArrayBox& yflux,FArrayBox& zflux),
                D_DECL(const int* ubc,const int* vbc,const int* wbc),
                D_DECL(const FArrayBox& U,const FArrayBox& V,const FArrayBox& W),
                D_DECL(int Ucomp,int Vcomp,int Wcomp),
                const FArrayBox& tforces, int Tcomp);

    // =============================================================
    // Data access functions
    // ==============================================================

    int useForcesInTrans() {return use_forces_in_trans;}
    
    // =============================================================
    // Advection functions follow
    // =============================================================

    void edge_states (const Box& grd, const Real* dx, Real dt, int velpred,
                      const FArrayBox& uedge, FArrayBox& stx,
                      const FArrayBox& vedge, FArrayBox& sty,
#if (BL_SPACEDIM == 3 )
                      const FArrayBox& wedge, FArrayBox& stz,
#endif
                      const FArrayBox& U, const FArrayBox& S, 
                      const FArrayBox& tforces, const FArrayBox& divu,
                      int fab_ind, int state_ind, const int *bc,
                      int use_conserv, AdvectionScheme which_scheme);

    void edge_states (const Box& grd, const Real* dx, Real dt, int velpred,
                      D_DECL(const FArrayBox &uedge, const FArrayBox &vedge, const FArrayBox &wedge),
                      D_DECL(int mCompX,       int mCompY,       int mCompZ),
                      D_DECL(FArrayBox &stx,   FArrayBox &sty,   FArrayBox &stz  ),
                      D_DECL(int eCompX,       int eCompY,       int eCompZ),
                      D_DECL(const FArrayBox& U, const FArrayBox &V, const FArrayBox &W),
                      D_DECL(int Ucomp,        int Vcomp,        int Wcomp),
                      const FArrayBox &S, int Scomp, const FArrayBox &tforces, int Tcomp, const FArrayBox& divu, int Dcomp,
                      int state_ind, const int *bc, int use_conserv, AdvectionScheme which_scheme);

    void edge_states_orig( const Box &grd, const Real *dx, Real dt, int velpred,
                           D_DECL(const FArrayBox &uedge, const FArrayBox &vedge, const FArrayBox &wedge),
                           D_DECL(int mCompX,       int mCompY,       int mCompZ),
                           D_DECL(FArrayBox &stx,   FArrayBox &sty,   FArrayBox &stz  ),
                           D_DECL(int eCompX,       int eCompY,       int eCompZ),
                           D_DECL(const FArrayBox& U, const FArrayBox &V, const FArrayBox &W),
                           D_DECL(int Ucomp,        int Vcomp,        int Wcomp),
                           const FArrayBox &S, int Scomp, const FArrayBox &tforces, int Tcomp,
                           int state_ind, const int *bc);

    void edge_states_fpu( const Box &grd, const Real *dx, Real dt,
                          D_DECL(const FArrayBox &uedge, const FArrayBox &vedge, const FArrayBox &wedge),
                          D_DECL(int mCompX,       int mCompY,       int mCompZ),
                          D_DECL(FArrayBox &stx,   FArrayBox &sty,   FArrayBox &stz  ),
                          D_DECL(int eCompX,       int eCompY,       int eCompZ),
                          const FArrayBox &S, int Scomp, const FArrayBox &tforces, int Tcomp, const FArrayBox& divu, int Dcomp,
                          int state_ind, const int *bc, int use_conserv);
    
    void edge_states_bds( const Box &grd, const Real *dx, Real dt,
                          D_DECL(const FArrayBox &uedge, const FArrayBox &vedge, const FArrayBox &wedge),
                          D_DECL(int mCompX,       int mCompY,       int mCompZ),
                          D_DECL(FArrayBox &stx,   FArrayBox &sty,   FArrayBox &stz  ),
                          D_DECL(int eCompX,       int eCompY,       int eCompZ),
                          const FArrayBox &S, int Scomp, const FArrayBox &tforces, int Tcomp, const FArrayBox& divu, int Dcomp,
                          int state_ind, const int *bc, int use_conserv);
    
    // compute divergent mac velocities
    void ComputeUmac( const Box &grd, const Real *dx, Real dt, 
                      FArrayBox &umac, const int *ubc, 
                      FArrayBox &vmac, const int *vbc, 
#if (BL_SPACEDIM == 3 )
                      FArrayBox &wmac, const int *wbc, 
#endif
                      FArrayBox &U, FArrayBox &tforces );
    
    // advect a state component
    void AdvectState( const Box &grd, const Real *dx, Real dt, 
                      FArrayBox &areax, FArrayBox &uedge, FArrayBox &xflux,  
                      FArrayBox &areay, FArrayBox &vedge, FArrayBox &yflux,  
#if (BL_SPACEDIM == 3 )                               
                      FArrayBox &areaz, FArrayBox &wedge, FArrayBox &zflux,
#endif
                      FArrayBox &U,
                      FArrayBox &S, 
                      FArrayBox &tforces, FArrayBox& divu,
                      int fab_ind,
                      FArrayBox &aofs,                  
                      int aofs_ind,
                      int iconserv, int state_ind, const int *bc,
                      AdvectionScheme scheme, FArrayBox &vol );

    // compute the advective derivative from fluxes
    void ComputeAofs (const Box& grd, 
                      const FArrayBox& areax, const FArrayBox& uedge, const FArrayBox& xflux,  
                      const FArrayBox& areay, const FArrayBox& vedge, const FArrayBox& yflux,  
#if (BL_SPACEDIM == 3 )
                      const FArrayBox& areaz, const FArrayBox& wedge, const FArrayBox& zflux,
#endif
                      const FArrayBox& vol,
                      FArrayBox& aofs, int aofs_ind, int iconserv) const;

    void ComputeAofs (const Box& grd, 
                      D_DECL(const FArrayBox& areax,
                             const FArrayBox& areay,
                             const FArrayBox& areaz),
                      D_DECL(int axcomp,int aycomp,int azcomp),
                      D_DECL(const FArrayBox& uedge,
                             const FArrayBox& vedge,
                             const FArrayBox& wedge),
                      D_DECL(int ucomp, int vcomp, int wcomp),
                      D_DECL(const FArrayBox& xflux,
                             const FArrayBox& yflux,
                             const FArrayBox& zflux),
                      D_DECL(int fxcomp,int fycomp,int fzcomp),
                      const FArrayBox& vol, int volcomp, 
                      FArrayBox& aofs,int acomp, int iconserv ) const;

    // sync advect a state component
    void SyncAdvect( const Box &grd, const Real *dx, Real dt, int level,
                     FArrayBox &areax, FArrayBox &uedge,
                     FArrayBox &ucorr, FArrayBox &xflux,
                     FArrayBox &areay, FArrayBox &vedge,
                     FArrayBox &vcorr, FArrayBox &yflux,
#if (BL_SPACEDIM == 3 )
                     FArrayBox &areaz, FArrayBox &wedge,
                     FArrayBox &wcorr, FArrayBox &zflux,
#endif
                     FArrayBox &U, FArrayBox &S,
                     FArrayBox &tforces, FArrayBox& divu, int ind,
                     FArrayBox &sync,                     int sync_ind,
                     int iconserv, int state_ind, const int *bc,
                     AdvectionScheme scheme, FArrayBox &vol );
    
    // compute the advective derivative of corrective fluxes for the mac sync
    void ComputeSyncAofs( const Box &grd,
                          FArrayBox &areax, FArrayBox &ucorr, FArrayBox &xflux,  
                          FArrayBox &areay, FArrayBox &vcorr, FArrayBox &yflux,  
#if (BL_SPACEDIM == 3 )
                          FArrayBox &areaz, FArrayBox &wcorr, FArrayBox &zflux,
#endif                     
                          FArrayBox &vol,
                          FArrayBox &sync,
                          int sync_ind, int iconserv );
    
    // correct a conservatively-advected scalar for under-over shoots
    void ConservativeScalMinMax( FArrayBox &Sold, FArrayBox &Snew, 
                                 int ind_old_s, int ind_old_rho,
                                 int ind_new_s, int ind_new_rho,
                                 const int *bc, const Box &grd );

    // correct a convectively-advected scalar for under-over shoots
    void ConvectiveScalMinMax( FArrayBox &Sold, FArrayBox &Snew, 
                               int ind_old, int ind_new,
                               const int *bc, const Box &grd );

    // =============================================================
    // Diagnostic functions follow
    // =============================================================
    
    // estimate the timestep from the cell-centered quantities
    Real estdt( FArrayBox &U, FArrayBox &tforces, FArrayBox &Rho,
                const Box &grd, const Real *dx, Real cfl, Real *u_max );

    // test the cell centered Courant number 
    Real test_u_rho( FArrayBox &U, FArrayBox &Rho, const Box &grd,
                     const Real *dx, const Real dt,
                     const Real *u_max );

    // test the mac edge velocity Courant numbers 
    Real test_umac_rho( FArrayBox &umac,
                        FArrayBox &vmac,
#if (BL_SPACEDIM == 3 )
                        FArrayBox &wmac,
#endif
                        FArrayBox &Rho,
                        const Box &grd, const Real *dx, const Real dt,
                        const Real *u_max );

    static bool are_any(const Array<AdvectionForm>& advectionType,
                        const AdvectionForm         testForm,
                        const int                   sComp,
                        const int                   nComp);

    static int how_many(const Array<AdvectionForm>& advectionType,
                        const AdvectionForm         testForm,
                        const int                   sComp,
                        const int                   nComp);

    static int hypgrow ();

    // =============================================================
    // Source term functions follow
    // =============================================================
    
    // compute the update rule, this is useful for 1st order RK
    // psi^n+1 = psi^n + dt*(tf^n)
    void Add_tf( const FArrayBox &Sold,
                 FArrayBox &Snew,    int start_ind, int num_comp,
                 const FArrayBox &tforces, int tf_ind,
                 const Box &grd,     Real dt ) const;

    // correct the 1st order RK to 2nd order via
    // psi^n+1 = psi^* + (dt/2)*(tf^* - tf^n)
    void Correct_tf( const FArrayBox &Sstar, FArrayBox &Snew,
                     int start_ind, int num_comp, 
                     const FArrayBox &tfstar, const FArrayBox &tfn,
                     int tf_ind,
                     const Box &grd,     Real dt ) const;

    // compute the update rule
    // psi^n+1 = psi^n - dt*aofs + dt*tforces
    void Add_aofs_tf( const FArrayBox &Sold,
                      FArrayBox &Snew,    int start_ind, int num_comp,
                      const FArrayBox &Aofs,    int aofs_ind,
                      const FArrayBox &tforces, int tf_ind,
                      const Box &grd,     Real dt ) const;

    // compute the update rule for velocities or momenta
    // psi^n+1 = psi^n - dt*aofs - dt*gp/rho + dt*tforces
    void Add_aofs_tf_gp( const FArrayBox &Uold, FArrayBox &Unew,
                         const FArrayBox &Aofs, const FArrayBox &tforces,
                         const FArrayBox &gp,   const FArrayBox &rho, 
                         const Box &grd,  Real dt) const;

    // compute total source term for velocities
    // tforces = (tforces + visc - gp)/rho
    void Sum_tf_gp_visc (FArrayBox&       tforces, int Tcomp,
                         const FArrayBox& visc,    int Vcomp,
                         const FArrayBox& gp,      int Gcomp,
                         const FArrayBox& rho,     int Rcomp) const;

    void Sum_tf_gp_visc( FArrayBox&       tforces,
                         const FArrayBox& visc, 
                         const FArrayBox& gp,
                         const FArrayBox& Rho) const;

    // iconserv==1   => tforces = tforces + visc - divU*S
    // iconserv==0   => tforces = (tforces+visc)/rho
    void Sum_tf_divu_visc (const FArrayBox& S,
                           FArrayBox& tforces,
                           int s_ind, int num_comp,
                           const FArrayBox& visc, int v_ind,
                           const FArrayBox& divu,
                           const FArrayBox& rho,
                           int iconserv) const;

    void Sum_tf_divu_visc (const FArrayBox& S, int s_ind,
                           FArrayBox& tforces, int t_ind,
                           int        num_comp,
                           const FArrayBox& visc,    int v_ind,
                           const FArrayBox& divu,    int d_ind,
                           const FArrayBox& rho,     int tr_ind,
                           int        iconserv) const;

    // compute total source term for velocities
    // tforces = (tforces - gp)/rho
    void Sum_tf_gp (FArrayBox& tforces, int Tcomp,
                    const FArrayBox& gp, int Gcomp,
                    const FArrayBox& rho, int Rcomp) const;

    // iconserv==1   => tforces = tforces - divU*S
    // iconserv==0   => tforces = (tforces)/rho
    void Sum_tf_divu (const FArrayBox& S,       int s_ind,
                      FArrayBox& tforces, int t_ind,
                      int        num_comp,
                      const FArrayBox& divu,    int d_ind,
                      const FArrayBox& rho,     int r_ind,
                      int        iconserv) const; 
    
    //================ for debugging purposes this is made public
    // internal FABS
    FArrayBox work; // work space for the Godunov Box
    FArrayBox uad;  // work space for advective transverse velocities
    FArrayBox vad;
#if (BL_SPACEDIM == 3 )    
    FArrayBox wad;
#endif

    // More workspace for modified advection scheme allowing for ppm option
    //  (pulled in to avoid f90 allocates)
    FArrayBox smp, dsvl, I;
    FArrayBox D_DECL(sedgex, sedgey, sedgez);

    // internal parameters
    int max_1d;    // largest dimension size allowed
    int scr_size;  // actual length of the scratch space
    
    // 1D arrays used in computing slopes and edges states
    Real *stxlo; 
    Real *stxhi; 
    Real *slxscr;
    Real *stylo;
    Real *styhi;
    Real *slyscr;
#if (BL_SPACEDIM == 3 )    
    Real *stzlo;
    Real *stzhi;
    Real *slzscr;
#endif
    
protected:   // protected work functions
    
    // manipulate 1d scratch space
    void ZeroScratch();
    void SetBogusScratch();
    void SetScratch( int max_size );
    void RemScratch();

    static void Initialize ();
    static void Finalize ();

protected:   // protected data members

    // sizing boxes
    Box work_bx;
    Box xflux_bx;
    Box yflux_bx;
#if (BL_SPACEDIM == 3 )    
    Box zflux_bx;
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

