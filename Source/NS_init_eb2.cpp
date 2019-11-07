
#include <NavierStokesBase.H>
#include <NAVIERSTOKESBASE_F.H>
using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include "AMReX_VisMF.H"
#include "AMReX_PlotFileUtil.H"



inline
bool NavierStokesBase::ebInitialized()
{
    return eb_initialized;
}

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
  
  const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());

  // These are the data sources
  volfrac = &(ebfactory.getVolFrac());
  bndrycent = &(ebfactory.getBndryCent());
  areafrac = ebfactory.getAreaFrac();
  facecent = ebfactory.getFaceCent();
  
  auto const& flags = ebfactory.getMultiEBCellFlagFab();

  for (MFIter mfi(*volfrac, false); mfi.isValid(); ++mfi)
  {
    BaseFab<int>& mfab = ebmask[mfi];
    const Box tbox = mfi.growntilebox();
    const FArrayBox& vfab = (*volfrac)[mfi];
    const EBCellFlagFab& flagfab = flags[mfi];
    
    FabType typ = flagfab.getType(tbox);
    int iLocal = mfi.LocalIndex();

    if (typ == FabType::regular) {
      mfab.setVal(1);
    }
    else if (typ == FabType::covered) {
      mfab.setVal(-1);
    }
    else if (typ == FabType::singlevalued) {
      int Ncut = 0;
      for (BoxIterator bit(tbox); bit.ok(); ++bit) {
        const EBCellFlag& flag = flagfab(bit(), 0);
    
        if (!(flag.isRegular() || flag.isCovered())) {
          Ncut++;
        }
      }
    
      for (BoxIterator bit(tbox); bit.ok(); ++bit) {
        const EBCellFlag& flag = flagfab(bit(), 0);
    
        if (!(flag.isRegular() || flag.isCovered())) {
          if (mfab.box().contains(bit())) mfab(bit()) = 0;
        } else {
          if (flag.isRegular()) {
            if (mfab.box().contains(bit())) mfab(bit()) = 1;
          } else if (flag.isCovered()) {
            if (mfab.box().contains(bit())) mfab(bit()) = -1;
          } else {
            if (mfab.box().contains(bit())) mfab(bit()) = 2;
          }
        }
      }
    }
    else {
      amrex::Print() << "unknown (or multivalued) fab type" << std::endl;
      amrex::Abort();
    }   
  }
}

void
NavierStokesBase::define_body_state()
{
  if (no_eb_in_domain) return;
  
  // Scan over data and find a point in the fluid to use to 
  // set computable values in cells outside the domain
  if (!body_state_set)
  {
    bool foundPt = false;
    const MultiFab& S = get_new_data(State_Type);
    BL_ASSERT(S.boxArray() == ebmask.boxArray());
    BL_ASSERT(S.DistributionMap() == ebmask.DistributionMap());
  
    body_state.resize(S.nComp(),0);
    for (MFIter mfi(S,false); mfi.isValid() && !foundPt; ++mfi)
    {
      const Box vbox = mfi.validbox();
      const BaseFab<int>& m = ebmask[mfi];
      const FArrayBox& fab = S[mfi];
      BL_ASSERT(m.box().contains(vbox));
  
      // TODO: Remove this dog and do this work in fortran 
      for (BoxIterator bit(vbox); bit.ok() && !foundPt; ++bit)
      {
        const IntVect& iv = bit();
        if (m(iv,0) == 1) {
          foundPt = true;
          for (int n=0; n<S.nComp(); ++n)
          {
            body_state[n] = fab(iv,n);
          }
        }
      }
    }
  
    // Find proc with lowest rank to find valid point, use that for all
    std::vector<int> found(ParallelDescriptor::NProcs(),0);
    found[ParallelDescriptor::MyProc()] = (int)foundPt;
    ParallelDescriptor::ReduceIntSum(&(found[0]),found.size());
    int body_rank = -1;
    for (int i=0; i<found.size(); ++i) {
      if (found[i]==1) {
        body_rank = i;
      }
    }
    BL_ASSERT(body_rank>=0);
    ParallelDescriptor::Bcast(&(body_state[0]),body_state.size(),body_rank);
    body_state_set = true;
  }
}

void
NavierStokesBase::set_body_state(MultiFab& S)
{
  if (no_eb_in_domain) return;

  if (!body_state_set)
  {
    define_body_state();
  }

  BL_ASSERT(S.nComp() == body_state.size());
  int nc = S.nComp();
  int covered_val = -1;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
  {
    const Box& vbox = mfi.validbox();
    fort_set_body_state(vbox.loVect(), vbox.hiVect(),
                      BL_TO_FORTRAN_ANYD(S[mfi]),
                      BL_TO_FORTRAN_ANYD(ebmask[mfi]),
                      &(body_state[0]),&nc,&covered_val);
  }
}



#endif
