#ifndef _GODUNOV_CLASS_H_
#define _GODUNOV_CLASS_H_ 

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

    // read parameters
    void read_params();
    
    // destruct the advection object
    ~Godunov();

    // set up the Godunov Box
    virtual void Setup( const BOX &grd, const REAL *dx, REAL dt, int set_fluxes, 
                        FARRAYBOX &xflux, int *ubc,
                        FARRAYBOX &yflux, int *vbc,
#if (BL_SPACEDIM == 3 )
                        FARRAYBOX &zflux, int *wbc,
#endif
                        FARRAYBOX &S, FARRAYBOX &Rho, 
                        const FARRAYBOX& tforces);

    // =============================================================
    // Data access functions
    // ==============================================================

    int useForcesInTrans() {return use_forces_in_trans;}
    
    // =============================================================
    // Advection functions follow
    // =============================================================

    // compute the edge states
    void edge_states( const BOX &grd, const REAL *dx, REAL dt, int velpred,
                      FARRAYBOX &uedge, FARRAYBOX &stx,
                      FARRAYBOX &vedge, FARRAYBOX &sty,
#if (BL_SPACEDIM == 3 )
                      FARRAYBOX &wedge, FARRAYBOX &stz,
#endif
                      FARRAYBOX &U,
                      FARRAYBOX &S, FARRAYBOX &tforces, int fab_ind,
                      int state_ind, int *bc);
    
    
    // compute divergent mac velocities
    void ComputeUmac( const BOX &grd, const REAL *dx, REAL dt, 
                      FARRAYBOX &umac, int *ubc, 
                      FARRAYBOX &vmac, int *vbc, 
#if (BL_SPACEDIM == 3 )
                      FARRAYBOX &wmac, int *wbc, 
#endif
                      FARRAYBOX &U, FARRAYBOX &tforces );
    
    // advect a state component
    virtual void AdvectState( const BOX &grd, const REAL *dx, REAL dt, 
                              FARRAYBOX &areax, FARRAYBOX &uedge, FARRAYBOX &xflux,  
                              FARRAYBOX &areay, FARRAYBOX &vedge, FARRAYBOX &yflux,  
#if (BL_SPACEDIM == 3 )                               
                              FARRAYBOX &areaz, FARRAYBOX &wedge, FARRAYBOX &zflux,
#endif
                              FARRAYBOX &U,
                              FARRAYBOX &S, FARRAYBOX &tforces, int fab_ind,
                              FARRAYBOX &aofs,                  int aofs_ind,
                              int iconserv, int state_ind, int *bc,
                              FARRAYBOX &vol );

    // compute the advective derivative from fluxes
    void ComputeAofs( const BOX &grd, 
                      FARRAYBOX &areax, FARRAYBOX &uedge, FARRAYBOX &xflux,  
                      FARRAYBOX &areay, FARRAYBOX &vedge, FARRAYBOX &yflux,  
#if (BL_SPACEDIM == 3 )
                      FARRAYBOX &areaz, FARRAYBOX &wedge, FARRAYBOX &zflux,
#endif
                      FARRAYBOX &vol,
                      FARRAYBOX &aofs,  int aofs_ind, int iconserv );

    // sync advect a state component
    virtual void SyncAdvect( const BOX &grd, const REAL *dx, REAL dt, int level,
                             FARRAYBOX &areax, FARRAYBOX &uedge,
                             FARRAYBOX &ucorr, FARRAYBOX &xflux,
                             FARRAYBOX &areay, FARRAYBOX &vedge,
                             FARRAYBOX &vcorr, FARRAYBOX &yflux,
#if (BL_SPACEDIM == 3 )
                             FARRAYBOX &areaz, FARRAYBOX &wedge,
                             FARRAYBOX &wcorr, FARRAYBOX &zflux,
#endif
                             FARRAYBOX &S, FARRAYBOX &tforces, int ind,
                             FARRAYBOX &sync,                  int sync_ind,
                             int iconserv, int state_ind, int *bc,
                             FARRAYBOX &vol );
    
    // compute the advective derivative of corrective fluxes for the mac sync
    void ComputeSyncAofs( const BOX &grd,
                          FARRAYBOX &areax, FARRAYBOX &ucorr, FARRAYBOX &xflux,  
                          FARRAYBOX &areay, FARRAYBOX &vcorr, FARRAYBOX &yflux,  
#if (BL_SPACEDIM == 3 )
                          FARRAYBOX &areaz, FARRAYBOX &wcorr, FARRAYBOX &zflux,
#endif                     
                          FARRAYBOX &vol,
                          FARRAYBOX &sync,
                          int sync_ind, int iconserv );
    
    // correct a scalar for under-over shoots
    void ScalMinMax( FARRAYBOX &Sold, FARRAYBOX &Snew, int ind, 
                     int *bc, const BOX &grd );

    // =============================================================
    // Diagnostic functions follow
    // =============================================================
    
    // estimate the timestep from the cell-centered quantities
    REAL estdt( FARRAYBOX &U, FARRAYBOX &tforces, FARRAYBOX &Rho,
                const BOX &grd, const REAL *dx, REAL cfl, REAL *u_max );

    // test the cell centered Courant number 
    REAL test_u_rho( FARRAYBOX &U, FARRAYBOX &Rho, const BOX &grd,
                     const REAL *dx, const REAL dt,
                     const REAL *u_max );

    // test the mac edge velocity Courant numbers 
    REAL test_umac_rho( FARRAYBOX &umac,
                        FARRAYBOX &vmac,
#if (BL_SPACEDIM == 3 )
                        FARRAYBOX &wmac,
#endif
                        FARRAYBOX &Rho,
                        const BOX &grd, const REAL *dx, const REAL dt,
                        const REAL *u_max );

    // =============================================================
    // Source term functions follow
    // =============================================================
    
    // compute the update rule, this is useful for 1st order RK
    // psi^n+1 = psi^n + dt*(tf^n)
    void Add_tf( FARRAYBOX &Sold,
                 FARRAYBOX &Snew,    int start_ind, int num_comp,
                 FARRAYBOX &tforces, int tf_ind,
                 const BOX &grd,     REAL dt );

    // correct the 1st order RK to 2nd order via
    // psi^n+1 = psi^* + (dt/2)*(tf^* - tf^n)
    void Correct_tf( FARRAYBOX &Sstar, FARRAYBOX &Snew,
                     int start_ind, int num_comp, 
                     FARRAYBOX &tfstar, FARRAYBOX &tfn,
                     int tf_ind,
                     const BOX &grd,     REAL dt );

    // compute the update rule
    // psi^n+1 = psi^n - dt*aofs + dt*tforces
    void Add_aofs_tf( FARRAYBOX &Sold,
                      FARRAYBOX &Snew,    int start_ind, int num_comp,
                      FARRAYBOX &Aofs,    int aofs_ind,
                      FARRAYBOX &tforces, int tf_ind,
                      const BOX &grd,     REAL dt );

    // compute the update rule for velocities
    // psi^n+1 = psi^n - dt*aofs - dt*gp/rho + dt*tforces
    void Add_aofs_tf_gp( FARRAYBOX &Uold, FARRAYBOX &Unew,
                         FARRAYBOX &Aofs, FARRAYBOX &tforces,
                         FARRAYBOX &gp,   FARRAYBOX &rho, 
                         const BOX &grd,  REAL dt );

    // compute total source term for velocities
    // tforces = (tforces + visc - gp)/rho
    void Sum_tf_gp_visc( FARRAYBOX &tforces, FARRAYBOX &visc, 
                         FARRAYBOX &gp,      FARRAYBOX &Rho );

    // iconserv==1   => tforces = tforces + visc - divU*S
    // iconserv==0   => tforces = (tforces+visc)/rho
    void Sum_tf_divu_visc( FARRAYBOX &S, FARRAYBOX &tforces,
                           int s_ind, int num_comp,
                           FARRAYBOX &visc,    int v_ind,
                           FARRAYBOX &divu,    FARRAYBOX &rho,
                           int iconserv );
    
    // compute total source term for velocities
    // tforces = (tforces - gp)/rho
    void Sum_tf_gp( FARRAYBOX &tforces, 
                    FARRAYBOX &gp,      FARRAYBOX &Rho );

    // iconserv==1   => tforces = tforces - divU*S
    // iconserv==0   => tforces = (tforces)/rho
    void Sum_tf_divu( FARRAYBOX &S, FARRAYBOX &tforces,
                      int s_ind, int num_comp,
                      FARRAYBOX &divu,    FARRAYBOX &rho,
                      int iconserv );

    //================ for debugging purposes this is made public
    // internal FABS
    FARRAYBOX work; // work space for the Godunov Box
    FARRAYBOX uad;  // work space for advective transverse velocities
    FARRAYBOX vad;
#if (BL_SPACEDIM == 3 )    
    FARRAYBOX wad;
#endif

    // internal parameters
    static REAL bogus_value;  // bogus value for debugging
    int max_1d;    // largest dimension size allowed
    int scr_size;  // actual length of the scratch space
    
    // 1D arrays used in computing slopes and edges states
    REAL *stxlo; 
    REAL *stxhi; 
    REAL *slxscr;
    REAL *stylo;
    REAL *styhi;
    REAL *slyscr;
#if (BL_SPACEDIM == 3 )    
    REAL *stzlo;
    REAL *stzhi;
    REAL *slzscr;
#endif
    
protected:   // protected work functions
    
    // manipulate 1d scratch space
    void ZeroScratch();
    void SetBogusScratch();
    void SetScratch( int max_size );
    void RemScratch();

protected:   // protected data members

    // sizing boxes
    BOX work_bx;
    BOX xflux_bx;
    BOX yflux_bx;
#if (BL_SPACEDIM == 3 )    
    BOX zflux_bx;
#endif

// flags

    static int printMinMax;
    static int use_forces_in_trans;

// algorithm control parameters

    static int slope_order;

};

#endif

