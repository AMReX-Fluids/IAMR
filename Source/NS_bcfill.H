#ifndef NS_bcfill_H_
#define NS_bcfill_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <NavierStokes.H>

using namespace amrex;

//
// Fill external Dirichlet boundary for State varibles.
// This simple version only makes use constant BC values set in the inputs file.
// Create a stuct, then use it to build a function that then gets passed to
// the StateDescriptor in NavierStokes::variableSetUp() NS_setup.cpp
//
struct stateFill
{
    amrex::GpuArray<amrex::GpuArray<amrex::Real, NavierStokes::NUM_STATE_MAX>, AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr stateFill ( amrex::GpuArray<amrex::GpuArray<amrex::Real, NavierStokes::NUM_STATE_MAX>,
			                  AMREX_SPACEDIM*2> const& a_bcv)
        : bcv(a_bcv) {}

// iv                  : Cell index
// dest, dcomp, numcomp: Fill numcomp components of dest starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// orig_comp           : component index for dcomp as in the desciptor set up in NS_setup.cpp 
  AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real /*time*/,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        using namespace amrex;

        // do something for external Dirichlet (BCType::ext_dir)

        const int i = iv[0];
        const int j = iv[1];
#if (AMREX_SPACEDIM == 3)
        const int k = iv[2];
#else
        const int k = 0;
#endif

        const Box& domain_box = geom.Domain();

        for (int nc = 0; nc < numcomp; ++nc)
	{
            const BCRec& bc = bcr[bcomp+nc];

            if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
            {
                dest(i,j,k,dcomp+nc) = bcv[Orientation(Direction::x,Orientation::low)][orig_comp+nc];
            }
            else if (bc.hi(0) == BCType::ext_dir and i > domain_box.bigEnd(0))
            {
                dest(i,j,k,dcomp+nc) = bcv[Orientation(Direction::x,Orientation::high)][orig_comp+nc];
            }

            if (bc.lo(1) == BCType::ext_dir and j < domain_box.smallEnd(1))
            {
                dest(i,j,k,dcomp+nc) = bcv[Orientation(Direction::y,Orientation::low)][orig_comp+nc];
            }
            else if (bc.hi(1) == BCType::ext_dir and j > domain_box.bigEnd(1))
            {
                dest(i,j,k,dcomp+nc) = bcv[Orientation(Direction::y,Orientation::high)][orig_comp+nc];
            }

#if (AMREX_SPACEDIM == 3)
            if (bc.lo(2) == BCType::ext_dir and k < domain_box.smallEnd(2))
            {
                dest(i,j,k,dcomp+nc) = bcv[Orientation(Direction::z,Orientation::low)][orig_comp+nc];
            }
            else if (bc.hi(2) == BCType::ext_dir and k > domain_box.bigEnd(2))
            {
                dest(i,j,k,dcomp+nc) = bcv[Orientation(Direction::z,Orientation::high)][orig_comp+nc];
            }
#endif
        }
    }
};

void state_fill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{
        GpuBndryFuncFab<stateFill> gpu_bndry_func(stateFill{NavierStokes::get_bc_values()});
        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);

}

//
// Fill external Dirichlet boundary for velocity
// Use this to create time and/or spatially dependent BCs
//
struct velFill
{
    int probtype;
  amrex::GpuArray<amrex::GpuArray<amrex::Real, NavierStokes::NUM_STATE_MAX>, AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr velFill (int a_probtype,
		       amrex::GpuArray<amrex::GpuArray<amrex::Real, NavierStokes::NUM_STATE_MAX>,
			                  AMREX_SPACEDIM*2> const& a_bcv)
      : probtype(a_probtype), bcv(a_bcv) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& vel,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real /*time*/,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        using namespace amrex;

	// do something for external Dirichlet (BCType::ext_dir)
		
        const int i = iv[0];
        const int j = iv[1];
#if (AMREX_SPACEDIM == 3)
        const int k = iv[2];
#else
        const int k = 0;
#endif

        const Box& domain_box = geom.Domain();

        for (int nc = 0; nc < numcomp; ++nc)
	{
            const BCRec& bc = bcr[bcomp+nc];

            if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
            {
		vel(i,j,k,dcomp+nc) = bcv[Orientation(Direction::x,Orientation::low)][orig_comp+nc];
		//
		// Do something interesting with velocity BC
		//
		// if (probtype ==  )
		// {
		//    // vel = some f(x,y,z,t) 
		// }
            }
            else if (bc.hi(0) == BCType::ext_dir and i > domain_box.bigEnd(0))
            {
                vel(i,j,k,dcomp+nc) = bcv[Orientation(Direction::x,Orientation::high)][orig_comp+nc];
            }

            if (bc.lo(1) == BCType::ext_dir and j < domain_box.smallEnd(1))
            {
                vel(i,j,k,dcomp+nc) = bcv[Orientation(Direction::y,Orientation::low)][orig_comp+nc];
            }
            else if (bc.hi(1) == BCType::ext_dir and j > domain_box.bigEnd(1))
            {
                vel(i,j,k,dcomp+nc) = bcv[Orientation(Direction::y,Orientation::high)][orig_comp+nc];
            }

#if (AMREX_SPACEDIM == 3)
            if (bc.lo(2) == BCType::ext_dir and k < domain_box.smallEnd(2))
            {
                vel(i,j,k,dcomp+nc) = bcv[Orientation(Direction::z,Orientation::low)][orig_comp+nc];
            }
            else if (bc.hi(2) == BCType::ext_dir and k > domain_box.bigEnd(2))
            {
                vel(i,j,k,dcomp+nc) = bcv[Orientation(Direction::z,Orientation::high)][orig_comp+nc];
            }
#endif
        }
    }
};

//
// Could add boundary Fill functions for density, tracer and temperature here.
// Make sure NavierStokes::variableSetUp() passes new Fill to StateDescriptor 
//

void vel_fill (Box const& bx, FArrayBox& data,
	       const int dcomp, const int numcomp,
	       Geometry const& geom, const Real time,
	       const Vector<BCRec>& bcr, const int bcomp,
	       const int scomp)
{

    GpuBndryFuncFab<velFill> gpu_bndry_func(velFill{NavierStokes::probtype,
	  NavierStokes::get_bc_values()});
    gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);

}

struct dummyFill
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    const int dcomp,
    const int numcomp,
    amrex::GeometryData const& geom,
    const amrex::Real /*time*/,
    const amrex::BCRec* bcr,
    const int bcomp,
    const int /*orig_comp*/) const
  {
    // Shouldn't actually ever use this, just need something computable.
    // Set to some ridiculous value so we know if it does get used.
    amrex::Real bogus_bc = 1.2345e40;

    const int i = iv[0];
    const int j = iv[1];
#if (AMREX_SPACEDIM == 3)
    const int k = iv[2];
#else
    const int k = 0;
#endif

    const Box& domain_box = geom.Domain();

    for (int nc = 0; nc < numcomp; ++nc)
    {
      const BCRec& bc = bcr[bcomp+nc];
      
      if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
      {
	dest(i,j,k,dcomp+nc) = bogus_bc;
	//
	// Do something interesting with destocity BC
	//
	// if (probtype ==  )
	// {
	//    // some f(x,y,z,t) to fill dest
	// }
      }
      else if (bc.hi(0) == BCType::ext_dir and i > domain_box.bigEnd(0))
      {
	dest(i,j,k,dcomp+nc) = bogus_bc;
      }
      
      if (bc.lo(1) == BCType::ext_dir and j < domain_box.smallEnd(1))
      {
	dest(i,j,k,dcomp+nc) = bogus_bc;
      }
      else if (bc.hi(1) == BCType::ext_dir and j > domain_box.bigEnd(1))
      {
	dest(i,j,k,dcomp+nc) = bogus_bc;
      }
#if (AMREX_SPACEDIM == 3)
      if (bc.lo(2) == BCType::ext_dir and k < domain_box.smallEnd(2))
      {
	dest(i,j,k,dcomp+nc) = bogus_bc;
      }
      else if (bc.hi(2) == BCType::ext_dir and k > domain_box.bigEnd(2))
      {
	dest(i,j,k,dcomp+nc) = bogus_bc;
      }
#endif
    }
  }
};

void dummy_fill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{

        GpuBndryFuncFab<dummyFill> gpu_bndry_func(dummyFill{});
        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);

}


// struct NodalFillExtDir
// {
//   AMREX_GPU_DEVICE
//   void operator()(
//     const amrex::IntVect& iv,
//     amrex::Array4<amrex::Real> const& dest,
//     const int dcomp,
//     const int numcomp,
//     amrex::GeometryData const& geom,
//     const amrex::Real time,
//     const amrex::BCRec* bcr,
//     const int bcomp,
//     const int orig_comp) const
//   {
//     // do something for external Dirichlet (BCType::ext_dir)
//   }
// };

void press_fill (Box const& /*bx*/, FArrayBox& /*data*/,
                 const int /*dcomp*/, const int /*numcomp*/,
                 Geometry const& /*geom*/, const Real /*time*/,
                 const Vector<BCRec>& /*bcr*/, const int /*bcomp*/,
                 const int /*scomp*/)
{
    amrex::Abort("press_fill: Need to write fill for external Dirichlet (BCType::ext_dir)");

    // GpuBndryFuncFab<NodalFillExtDir> gpu_bndry_func(NodalFillExtDir{});
    // gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);

}

#endif
