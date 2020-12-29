
#ifndef NS_DERIVE_CPP_H
#define NS_DERIVE_CPP_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>


    void der_vel_avg (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                      const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                      amrex::Real time, const int* bcrec, int level);


#endif
