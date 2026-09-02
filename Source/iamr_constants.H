#ifndef IAMR_CONSTANTS_H_
#define IAMR_CONSTANTS_H_

#include <AMReX_REAL.H>

#include <numbers>

// Value to use in covered cells
static constexpr amrex::Real COVERED_VAL = 1.0e40;

static constexpr amrex::Real Pi    = std::numbers::pi_v<amrex::Real>;
static constexpr amrex::Real TwoPi = amrex::Real(2.0) * std::numbers::pi_v<amrex::Real>;

#endif
