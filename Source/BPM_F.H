#ifndef _BPM_F_H_
#define _BPM_F_H_

/*
 * C prototypes and structures to match the C callable Fortran routines in bpm.f90.
 */

extern "C" {

  typedef struct bpm_parameters {
    int    type;
    double diameter;
  } bpm_parameters;

  void bpm_compute_forces_2d(double const *x1, double const *v1, double *f1,
                             double const *x2, double const *v2, double *f2,
                             double sconstant, double slength);

  void bpm_spread_forces_2d(double const *x1, double const *f1, bpm_parameters *params,
                            int const *lo, int const *hi, int ng, double const *xlo, double const *xhi, double dx,
                            double *f, double *divf);
}


#endif
