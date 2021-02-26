
#ifndef NS_DERIVE_CPP_H
#define NS_DERIVE_CPP_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>


namespace derive_functions {

  void der_vel_avg (const amrex::Box& bx,
		    amrex::FArrayBox& derfab, int dcomp, int ncomp,
		    const amrex::FArrayBox& datafab,
		    const amrex::Geometry& geomdata,
		    amrex::Real time, const int* bcrec, int level);

  //
  //  Compute cell-centered pressure as average of the 
  //  surrounding nodal values 
  //
  void deravgpres (const amrex::Box& bx,
		   amrex::FArrayBox& derfab, int dcomp, int ncomp,
		   const amrex::FArrayBox& datafab,
		   const amrex::Geometry& geomdata,
		   amrex::Real time, const int* bcrec, int level);

  //
  // Compute the magnitude of the vorticity
  //
  void dermgvort (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
		  const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
		  amrex::Real time, const int* bcrec, int level);
  //
  // Compute kinetic energy
  //
  void derkeng (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
		const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
		amrex::Real time, const int* bcrec, int level);
  //
  // Do nothing
  //
  void dernull (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
		const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
		amrex::Real time, const int* bcrec, int level);

}
#endif
