#include <NavierStokesBase.H>
#include "NS_derive.H"
#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#endif

using namespace amrex;

namespace derive_functions
{
  void der_vel_avg (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
		    const FArrayBox& datfab, const Geometry& /*geomdata*/,
		    Real time, const int* /*bcrec*/, int level)

  {
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= AMREX_SPACEDIM*2);
    AMREX_ASSERT(ncomp == AMREX_SPACEDIM*2);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
    amrex::Real inv_time;
    amrex::Real inv_time_fluct;

    if (NavierStokesBase::time_avg[level] == 0){
      inv_time = 1.0;
    }else{
      inv_time = 1.0 / NavierStokesBase::time_avg[level];
    }

    if (NavierStokesBase::time_avg_fluct[level] == 0){
      inv_time_fluct = 1.0;
    }else{
      inv_time_fluct = 1.0 / NavierStokesBase::time_avg_fluct[level];
    }

    amrex::ParallelFor(bx, AMREX_SPACEDIM, [inv_time,inv_time_fluct,der,in_dat]
		       AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
		       {
			 der(i,j,k,n) = in_dat(i,j,k,n) * inv_time;
			 der(i,j,k,n+AMREX_SPACEDIM) = sqrt(in_dat(i,j,k,n+AMREX_SPACEDIM) * inv_time_fluct);
		       });
  }

  //
  //  Compute cell-centered pressure as average of the 
  //  surrounding nodal values 
  //
  void deravgpres (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
		   const FArrayBox& datfab, const Geometry& /*geomdata*/,
		   Real /*time*/, const int* /*bcrec*/, int /*level*/)

  {
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(Box(datfab.box()).enclosedCells().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= 1);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);
#if (AMREX_SPACEDIM == 2 )
    Real factor = 0.25;
#elif (AMREX_SPACEDIM == 3 )
    Real factor = 0.125;
#endif

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      der(i,j,k) =  factor * (  in_dat(i+1,j,k)     + in_dat(i,j,k)
				+ in_dat(i+1,j+1,k)   + in_dat(i,j+1,k)
#if (AMREX_SPACEDIM == 3 )
				+ in_dat(i+1,j,k+1)   + in_dat(i,j,k+1)
				+ in_dat(i+1,j+1,k+1) + in_dat(i,j+1,k+1)
#endif
				);
    });
  }
  
  //
  //  Compute magnitude of vorticity
  //
  void dermgvort (const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
		  const FArrayBox& datfab, const Geometry& geomdata,
		  Real /*time*/, const int* /*bcrec*/, int /*level*/)

  {
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));

    AMREX_D_TERM(const amrex::Real idx = geomdata.InvCellSize(0);,
                 const amrex::Real idy = geomdata.InvCellSize(1);,
                 const amrex::Real idz = geomdata.InvCellSize(2););

    amrex::Array4<amrex::Real const> const& dat_arr = datfab.const_array();
    amrex::Array4<amrex::Real>       const&vort_arr = derfab.array();

#ifdef AMREX_USE_EB
    const EBFArrayBox& ebfab = static_cast<EBFArrayBox const&>(datfab);
    const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();
    auto typ = flags.getType(bx);
    if (typ == FabType::covered)
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
			   {
			     vort_arr(i,j,k) = 0.0;
			   });
    } else if (typ == FabType::singlevalued)
    {
	const auto& flag_fab = flags.const_array();
	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
	  constexpr amrex::Real c0 = -1.5;
	  constexpr amrex::Real c1 = 2.0;
	  constexpr amrex::Real c2 = -0.5;
	  if (flag_fab(i,j,k).isCovered()) {
	    vort_arr(i,j,k) = 0.0;
	  } else {
	    amrex::Real vx = 0.0;
	    amrex::Real uy = 0.0;
#if ( AMREX_SPACEDIM == 2 )
	    // Need to check if there are covered cells in neighbours --
	    // -- if so, use one-sided difference computation (but still quadratic)
	    if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
	      vx = - (c0 * dat_arr(i  ,j,k,1)
		      + c1 * dat_arr(i-1,j,k,1)
		      + c2 * dat_arr(i-2,j,k,1)) * idx;
	    } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
	      vx = (c0 * dat_arr(i  ,j,k,1)
		    + c1 * dat_arr(i+1,j,k,1)
		    + c2 * dat_arr(i+2,j,k,1)) * idx;
	    } else {
	      vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
	    }
	    // Do the same in y-direction
	    if (!flag_fab(i,j,k).isConnected( 0,1,0)) {
	      uy = - (c0 * dat_arr(i,j  ,k,0)
		      + c1 * dat_arr(i,j-1,k,0)
		      + c2 * dat_arr(i,j-2,k,0)) * idy;
	    } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
	      uy = (c0 * dat_arr(i,j  ,k,0)
		    + c1 * dat_arr(i,j+1,k,0)
		    + c2 * dat_arr(i,j+2,k,0)) * idy;
	    } else {
	      uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
	    }
	    vort_arr(i,j,k) = vx-uy;
#elif ( AMREX_SPACEDIM == 3 )
	    amrex::Real wx = 0.0;
	    amrex::Real wy = 0.0;
	    amrex::Real uz = 0.0;
	    amrex::Real vz = 0.0;
	    // Need to check if there are covered cells in neighbours --
	    // -- if so, use one-sided difference computation (but still quadratic)
	    if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
	      // Covered cell to the right, go fish left
	      vx = - (c0 * dat_arr(i  ,j,k,1)
		      + c1 * dat_arr(i-1,j,k,1)
		      + c2 * dat_arr(i-2,j,k,1)) * idx;
	      wx = - (c0 * dat_arr(i  ,j,k,2)
		      + c1 * dat_arr(i-1,j,k,2)
		      + c2 * dat_arr(i-2,j,k,2)) * idx;
	    } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
	      // Covered cell to the left, go fish right
	      vx = (c0 * dat_arr(i  ,j,k,1)
		    + c1 * dat_arr(i+1,j,k,1)
		    + c2 * dat_arr(i+2,j,k,1)) * idx;
	      wx = (c0 * dat_arr(i  ,j,k,2)
		    + c1 * dat_arr(i+1,j,k,2)
		    + c2 * dat_arr(i+2,j,k,2)) * idx;
	    } else {
	      // No covered cells right or left, use standard stencil
	      vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
	      wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;
	    }
	    // Do the same in y-direction
	    if (!flag_fab(i,j,k).isConnected(0, 1,0)) {
	      uy = - (c0 * dat_arr(i,j  ,k,0)
		      + c1 * dat_arr(i,j-1,k,0)
		      + c2 * dat_arr(i,j-2,k,0)) * idy;
	      wy = - (c0 * dat_arr(i,j  ,k,2)
		      + c1 * dat_arr(i,j-1,k,2)
		      + c2 * dat_arr(i,j-2,k,2)) * idy;
	    } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
	      uy = (c0 * dat_arr(i,j  ,k,0)
		    + c1 * dat_arr(i,j+1,k,0)
		    + c2 * dat_arr(i,j+2,k,0)) * idy;
	      wy = (c0 * dat_arr(i,j  ,k,2)
		    + c1 * dat_arr(i,j+1,k,2)
		    + c2 * dat_arr(i,j+2,k,2)) * idy;
	    } else {
	      uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
	      wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;
	    }
	    // Do the same in z-direction
	    if (!flag_fab(i,j,k).isConnected(0,0, 1)) {
	      uz = - (c0 * dat_arr(i,j,k  ,0)
		      + c1 * dat_arr(i,j,k-1,0)
		      + c2 * dat_arr(i,j,k-2,0)) * idz;
	      vz = - (c0 * dat_arr(i,j,k  ,1)
		      + c1 * dat_arr(i,j,k-1,1)
		      + c2 * dat_arr(i,j,k-2,1)) * idz;
	    } else if (!flag_fab(i,j,k).isConnected(0,0,-1)) {
	      uz = (c0 * dat_arr(i,j,k  ,0)
		    + c1 * dat_arr(i,j,k+1,0)
		    + c2 * dat_arr(i,j,k+2,0)) * idz;
	      vz = (c0 * dat_arr(i,j,k  ,1)
		    + c1 * dat_arr(i,j,k+1,1)
		    + c2 * dat_arr(i,j,k+2,1)) * idz;
	    } else {
	      uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
	      vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;
	    }
	    vort_arr(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
	  }
	});
    } else
#endif
      {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if ( AMREX_SPACEDIM == 2 )
	  amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
	  amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
	  vort_arr(i,j,k) = amrex::Math::abs(vx-uy);
	  
#elif ( AMREX_SPACEDIM == 3 )
	  amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;

	  amrex::Real wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;
	  
	  amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
	  amrex::Real wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;
	  
	  amrex::Real uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
	  amrex::Real vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;
	  
	  vort_arr(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
	});
      }
  }

  //
  // Compute kinetic energy
  //
  void derkeng (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
		const FArrayBox& datfab, const Geometry& /*geomdata*/,
		Real /*time*/, const int* /*bcrec*/, int /*level*/)

  {
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(Box(datfab.box()).contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= 1);
    AMREX_ASSERT(ncomp == 1);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);

    amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const Real rho = in_dat(i,j,k,0);
      const Real vx  = in_dat(i,j,k,1);
      const Real vy  = in_dat(i,j,k,2);
#if (AMREX_SPACEDIM ==3)
      const Real vz  = in_dat(i,j,k,3);
#endif
      
      der(i,j,k) =  0.5 * rho * ( vx*vx + vy*vy 
#if (AMREX_SPACEDIM == 3 )
				  + vz*vz
#endif
				  );
    });
  }

  //
  // Null function
  //
  void dernull (const Box& /*bx*/,
		FArrayBox& /*derfab*/, int /*dcomp*/, int /*ncomp*/,
		const FArrayBox& /*datfab*/, const Geometry& /*geomdata*/,
		Real /*time*/, const int* /*bcrec*/, int /*level*/)

  {
    //
    // Do nothing. 
    //
  }

};
