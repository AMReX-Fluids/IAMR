#ifndef _GODUNOV_CLASS_H_
#define _GODUNOV_CLASS_H_ 

//
// $Id: Godunov.H,v 1.3 1997-10-01 01:03:07 car Exp $
//

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
                        FArrayBox &xflux, int *ubc,
                        FArrayBox &yflux, int *vbc,
#if (BL_SPACEDIM == 3 )
                        FArrayBox &zflux, int *wbc,
#endif
                        FArrayBox &S, FArrayBox &Rho, 
                        const FArrayBox& tforces);

    // =============================================================
    // Data access functions
    // ==============================================================

    int useForcesInTrans() {return use_forces_in_trans;}
    
    // =============================================================
    // Advection functions follow
    // =============================================================

    // compute the edge states
    void edge_states( const BOX &grd, const REAL *dx, REAL dt, int velpred,
                      FArrayBox &uedge, FArrayBox &stx,
                      FArrayBox &vedge, FArrayBox &sty,
#if (BL_SPACEDIM == 3 )
                      FArrayBox &wedge, FArrayBox &stz,
#endif
                      FArrayBox &U,
                      FArrayBox &S, FArrayBox &tforces, int fab_ind,
                      int state_ind, int *bc);
    
    
    // compute divergent mac velocities
    void ComputeUmac( const BOX &grd, const REAL *dx, REAL dt, 
                      FArrayBox &umac, int *ubc, 
                      FArrayBox &vmac, int *vbc, 
#if (BL_SPACEDIM == 3 )
                      FArrayBox &wmac, int *wbc, 
#endif
                      FArrayBox &U, FArrayBox &tforces );
    
    // advect a state component
    virtual void AdvectState( const BOX &grd, const REAL *dx, REAL dt, 
                              FArrayBox &areax, FArrayBox &uedge, FArrayBox &xflux,  
                              FArrayBox &areay, FArrayBox &vedge, FArrayBox &yflux,  
#if (BL_SPACEDIM == 3 )                               
                              FArrayBox &areaz, FArrayBox &wedge, FArrayBox &zflux,
#endif
                              FArrayBox &U,
                              FArrayBox &S, FArrayBox &tforces, int fab_ind,
                              FArrayBox &aofs,                  int aofs_ind,
                              int iconserv, int state_ind, int *bc,
                              FArrayBox &vol );

    // compute the advective derivative from fluxes
    void ComputeAofs( const BOX &grd, 
                      FArrayBox &areax, FArrayBox &uedge, FArrayBox &xflux,  
                      FArrayBox &areay, FArrayBox &vedge, FArrayBox &yflux,  
#if (BL_SPACEDIM == 3 )
                      FArrayBox &areaz, FArrayBox &wedge, FArrayBox &zflux,
#endif
                      FArrayBox &vol,
                      FArrayBox &aofs,  int aofs_ind, int iconserv );

    // sync advect a state component
    virtual void SyncAdvect( const BOX &grd, const REAL *dx, REAL dt, int level,
                             FArrayBox &areax, FArrayBox &uedge,
                             FArrayBox &ucorr, FArrayBox &xflux,
                             FArrayBox &areay, FArrayBox &vedge,
                             FArrayBox &vcorr, FArrayBox &yflux,
#if (BL_SPACEDIM == 3 )
                             FArrayBox &areaz, FArrayBox &wedge,
                             FArrayBox &wcorr, FArrayBox &zflux,
#endif
                             FArrayBox &S, FArrayBox &tforces, int ind,
                             FArrayBox &sync,                  int sync_ind,
                             int iconserv, int state_ind, int *bc,
                             FArrayBox &vol );
    
    // compute the advective derivative of corrective fluxes for the mac sync
    void ComputeSyncAofs( const BOX &grd,
                          FArrayBox &areax, FArrayBox &ucorr, FArrayBox &xflux,  
                          FArrayBox &areay, FArrayBox &vcorr, FArrayBox &yflux,  
#if (BL_SPACEDIM == 3 )
                          FArrayBox &areaz, FArrayBox &wcorr, FArrayBox &zflux,
#endif                     
                          FArrayBox &vol,
                          FArrayBox &sync,
                          int sync_ind, int iconserv );
    
    // correct a scalar for under-over shoots
    void ScalMinMax( FArrayBox &Sold, FArrayBox &Snew, int ind, 
                     int *bc, const BOX &grd );

    // =============================================================
    // Diagnostic functions follow
    // =============================================================
    
    // estimate the timestep from the cell-centered quantities
    REAL estdt( FArrayBox &U, FArrayBox &tforces, FArrayBox &Rho,
                const BOX &grd, const REAL *dx, REAL cfl, REAL *u_max );

    // test the cell centered Courant number 
    REAL test_u_rho( FArrayBox &U, FArrayBox &Rho, const BOX &grd,
                     const REAL *dx, const REAL dt,
                     const REAL *u_max );

    // test the mac edge velocity Courant numbers 
    REAL test_umac_rho( FArrayBox &umac,
                        FArrayBox &vmac,
#if (BL_SPACEDIM == 3 )
                        FArrayBox &wmac,
#endif
                        FArrayBox &Rho,
                        const BOX &grd, const REAL *dx, const REAL dt,
                        const REAL *u_max );

    // =============================================================
    // Source term functions follow
    // =============================================================
    
    // compute the update rule, this is useful for 1st order RK
    // psi^n+1 = psi^n + dt*(tf^n)
    void Add_tf( FArrayBox &Sold,
                 FArrayBox &Snew,    int start_ind, int num_comp,
                 FArrayBox &tforces, int tf_ind,
                 const BOX &grd,     REAL dt );

    // correct the 1st order RK to 2nd order via
    // psi^n+1 = psi^* + (dt/2)*(tf^* - tf^n)
    void Correct_tf( FArrayBox &Sstar, FArrayBox &Snew,
                     int start_ind, int num_comp, 
                     FArrayBox &tfstar, FArrayBox &tfn,
                     int tf_ind,
                     const BOX &grd,     REAL dt );

    // compute the update rule
    // psi^n+1 = psi^n - dt*aofs + dt*tforces
    void Add_aofs_tf( FArrayBox &Sold,
                      FArrayBox &Snew,    int start_ind, int num_comp,
                      FArrayBox &Aofs,    int aofs_ind,
                      FArrayBox &tforces, int tf_ind,
                      const BOX &grd,     REAL dt );

    // compute the update rule for velocities
    // psi^n+1 = psi^n - dt*aofs - dt*gp/rho + dt*tforces
    void Add_aofs_tf_gp( FArrayBox &Uold, FArrayBox &Unew,
                         FArrayBox &Aofs, FArrayBox &tforces,
                         FArrayBox &gp,   FArrayBox &rho, 
                         const BOX &grd,  REAL dt );

    // compute total source term for velocities
    // tforces = (tforces + visc - gp)/rho
    void Sum_tf_gp_visc( FArrayBox &tforces, FArrayBox &visc, 
                         FArrayBox &gp,      FArrayBox &Rho );

    // iconserv==1   => tforces = tforces + visc - divU*S
    // iconserv==0   => tforces = (tforces+visc)/rho
    void Sum_tf_divu_visc( FArrayBox &S, FArrayBox &tforces,
                           int s_ind, int num_comp,
                           FArrayBox &visc,    int v_ind,
                           FArrayBox &divu,    FArrayBox &rho,
                           int iconserv );
    
    // compute total source term for velocities
    // tforces = (tforces - gp)/rho
    void Sum_tf_gp( FArrayBox &tforces, 
                    FArrayBox &gp,      FArrayBox &Rho );

    // iconserv==1   => tforces = tforces - divU*S
    // iconserv==0   => tforces = (tforces)/rho
    void Sum_tf_divu( FArrayBox &S, FArrayBox &tforces,
                      int s_ind, int num_comp,
                      FArrayBox &divu,    FArrayBox &rho,
                      int iconserv );

    //================ for debugging purposes this is made public
    // internal FABS
    FArrayBox work; // work space for the Godunov Box
    FArrayBox uad;  // work space for advective transverse velocities
    FArrayBox vad;
#if (BL_SPACEDIM == 3 )    
    FArrayBox wad;
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

