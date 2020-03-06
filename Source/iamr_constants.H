#ifndef IAMR_CONSTANTS_H_
#define IAMR_CONSTANTS_H_

#include <AMReX_REAL.H>
#include <cmath>

// What we consider a ~zero velocity when upwinding
static constexpr amrex::Real small_vel = 1.e-10;

#endif
