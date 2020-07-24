#ifndef IAMR_GODUNOV_CONSTANTS_H_
#define IAMR_GODUNOV_CONSTANTS_H_

#include <AMReX_REAL.H>
#include <cmath>


#define XVEL 0
#define YVEL 1
#define ZVEL 2

// What we consider a ~zero velocity when upwinding
static constexpr amrex::Real small_vel = 1.e-8;


#endif
