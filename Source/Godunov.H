
#ifndef _GODUNOV_CLASS_H_
#define _GODUNOV_CLASS_H_

#define MICHELE_DEBUG 0

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

    // flags

    static int verbose;
    static int use_forces_in_trans;

// algorithm control parameters

    static int  slope_order;
    static int  ppm_type;
    static int  corner_couple;

    // =============================================================
    // Setup functions follow
    // =============================================================

    // construct the advection object
    Godunov(int max_size = 512);


    // destruct the advection object
    ~Godunov();


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

    static int hypgrow ();


protected:   // protected work functions

    static void Initialize ();
    static void Finalize ();


};

#endif
