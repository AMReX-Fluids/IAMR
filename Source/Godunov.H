
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

//     static int verbose;
//     static int use_forces_in_trans;

// // algorithm control parameters

//     static int  slope_order;
//     static int  ppm_type;
//     static int  corner_couple;

    // =============================================================
    // Setup functions follow
    // =============================================================

    // construct the advection object
    Godunov(int max_size = 512);


    // destruct the advection object
    ~Godunov();

    // =============================================================
    // Diagnostic functions follow
    // =============================================================

    static int hypgrow ();


protected:   // protected work functions

    static void Initialize ();
    static void Finalize ();

    static int verbose;
    static int use_forces_in_trans;

// algorithm control parameters

    static int  slope_order;
    static int  ppm_type;
    static int  corner_couple;


};

#endif
