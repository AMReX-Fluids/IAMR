




#include <NavierStokesBase.H>
#include <AMReX_VisMF.H>



#include <AMReX_MLMG.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLEBTensorOp.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB_utils.H>
#else
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLTensorOp.H>
#endif


using namespace amrex;



void
NavierStokesBase::calc_mut_LES(const Real time)

{

  if (ParallelDescriptor::IOProcessor() && getLESVerbose) {
    amrex::Print() << "\n DEBUG WE ARE IN THE NEW ROUTINE FOR LES \n\n";
  }
  
  auto whichTime = which_time(State_Type,time);
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  MultiFab& Sstate = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);
  
  int nGrow = 1;
  FillPatchIterator fpi(*this,Sstate,nGrow,time,State_Type,Xvel,AMREX_SPACEDIM);
  MultiFab& Uvel=fpi.get_mf();
  


      int agglomeration = 1;
      int consolidation = 1;
      LPInfo info;
      info.setAgglomeration(agglomeration);
      info.setConsolidation(consolidation);
      //fixme??
      info.setMetricTerm(false);
      info.setMaxCoarseningLevel(100);

#ifdef AMREX_USE_EB
      const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory());
      MLEBTensorOp tensorop({geom}, {grids}, {dmap}, info, {ebf});

      std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM>areafrac = ebf->getAreaFrac();
#else
      MLTensorOp tensorop({geom}, {grids}, {dmap}, info);
#endif


      int tensor_max_order    = 2;
      tensorop.setMaxOrder(tensor_max_order);




      // create right container
      Array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc[AMREX_SPACEDIM];
      Array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc[AMREX_SPACEDIM];
      // fill it
      for (int i=0; i<AMREX_SPACEDIM; i++){
        LES_setDomainBC(mlmg_lobc[i], mlmg_hibc[i], Xvel+i);
      }    
        
      // pass to op
      tensorop.setDomainBC({AMREX_D_DECL(mlmg_lobc[0],mlmg_lobc[1],mlmg_lobc[2])},
                           {AMREX_D_DECL(mlmg_hibc[0],mlmg_hibc[1],mlmg_hibc[2])});  
  
     
      
      // set up level BCs
      // WARNING THIS PART BELOW HAS NOT BEEN TESTED
      //{
      //  MultiFab crsedata;
      //  int ng = 0;
      //  const int soln_ng = 1;
      //
      //  if (level > 0) {
      //    //auto& crse_ns = *(coarser->navier_stokes);
      //    NavierStokesBase& crse_ns  = getLevel(level-1);
      //    crsedata.define(crse_ns.boxArray(), crse_ns.DistributionMap(), AMREX_SPACEDIM,
      //                    ng, MFInfo(),crse_ns.Factory());
      //    AmrLevel::FillPatch(crse_ns, crsedata, ng, time, State_Type, Xvel,
      //                        AMREX_SPACEDIM);
      //    tensorop.setCoarseFineBC(&crsedata, crse_ratio[0]);
      //  }
      //
      //  AmrLevel::FillPatch(*this,Uvel,soln_ng,time,State_Type,Xvel,AMREX_SPACEDIM);
      //
        tensorop.setLevelBC(0, &Uvel);
      //}
      
      
        tensorop.setScalars(0., 1.0);
        tensorop.setACoeffs(0, 0.0);
        tensorop.setShearViscosity(0, [1. ,1.]);
        tensorop.setBulkViscosity(0, 0.);
#ifdef AMREX_USE_EB
        tensorop.setEBShearViscosity(0, 1.);
        tensorop.setEBBulkViscosity(0, 0.);
#endif




 FluxBoxes fb(this,AMREX_SPACEDIM);
 MultiFab** tensorflux = fb.get();
std::array<MultiFab*,AMREX_SPACEDIM> fp{AMREX_D_DECL(tensorflux[0], tensorflux[1], tensorflux[2])};

  // OK READY TO CREATE THE ROUTINE TO COMPUTE GRADIENT INSTEAD OF CALLING COMPFLUX 
 tensorop.compFlux(0,{fp},{Uvel},MLLinOp::Location::FaceCenter);
 
 //VisMF::Write(Uvel,"test");
 
  VisMF::Write(*fp[0],"test");

  
  if (LES_model == "Smagorinsky") {
    amrex::Print() << "\n DEBUG WE DO SMAGORNISKY \n\n";
  }
  else
  {
    amrex::Abort("\n DEBUG DONT KNOW THIS LES MODEL \n\n");
  }
  
}




void
NavierStokesBase::LES_setDomainBC (std::array<LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
                        std::array<LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc,
                        int src_comp)
{
    const BCRec& bc = get_desc_lst()[State_Type].getBC(src_comp);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (parent->Geom(0).isPeriodic(idim))
        {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else
        {
            int pbc = bc.lo(idim);
            if (pbc == EXT_DIR)
            {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      ||
                     pbc == REFLECT_EVEN)
            {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }
            else if (pbc == REFLECT_ODD)
            {
                mlmg_lobc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                mlmg_lobc[idim] = LinOpBCType::bogus;
            }

            pbc = bc.hi(idim);
            if (pbc == EXT_DIR)
            {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      ||
                     pbc == REFLECT_EVEN)
            {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
            else if (pbc == REFLECT_ODD)
            {
                mlmg_hibc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                mlmg_hibc[idim] = LinOpBCType::bogus;
            }
        }
    }
}
