#include <NavierStokesBase.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include "AMReX_VisMF.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

// called in main before Amr->init(start,stop)
void
initialize_EB2 (const Geometry& geom, const int required_coarsening_level,
		const int max_coarsening_level)
{
    // read in EB parameters
    ParmParse ppeb2("eb2");
    std::string geom_type;
    ppeb2.get("geom_type", geom_type);

    Print()<<"geom "<<geom<<"\n\n";
    //Build geometry -- WIP!!!
    // seems like this is a call that can create the implicit function,
    // GeometryShop, and Geometry all at "once"... geom might already be init
    //AMReX_EB2::Build (const Geometry& geom, int required_coarsening_level,
    //   int max_coarsening_level, int ngrow)
    //fixme -- not sure that 4 is the right number here, just taken from CNS
    EB2::Build(geom, required_coarsening_level, max_coarsening_level, 4);
}

void
NavierStokesBase::init_eb (const Geometry& level_geom, const BoxArray& ba, const DistributionMapping& dm)
{
  // Build the geometry information; this is done for each new set of grids
  initialize_eb2_structs();

}

void
NavierStokesBase::initialize_eb2_structs() {

  amrex::Print() << "Initializing EB2 structs" << std::endl;

  // NOTE: THIS NEEDS TO BE REPLACED WITH A FLAGFAB 
  
  // n.b., could set this to 1 if geometry is all_regular as an optimization
  no_eb_in_domain = 0;

  //  1->regular, 0->irregular, -1->covered, 2->outside
  ebmask.define(grids, dmap,  1, 0);
  
  //static_assert(std::is_standard_layout<EBBndryGeom>::value,
  //              "EBBndryGeom is not standard layout");

  //const amrex::MultiFab* volfrac;
  //const amrex::MultiCutFab* bndrycent;
  //std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> eb2areafrac;
  //std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> facecent;

  const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());

  // These are the data sources
  volfrac = &(ebfactory.getVolFrac());
  bndrycent = &(ebfactory.getBndryCent());
  areafrac = ebfactory.getAreaFrac();
  facecent = ebfactory.getFaceCent();
  
}

#endif
