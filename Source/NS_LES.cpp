




#include <NavierStokesBase.H>


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

  MultiFab& U = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);
  
  int nGrow = 1;
  FillPatchIterator fpi(*this,U,nGrow,time,State_Type,Xvel,AMREX_SPACEDIM);
  
  


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
      for (int i=0; i<AMREX_SPACEDIM; i++)
        Diffusion::setDomainBC(mlmg_lobc[i], mlmg_hibc[i], Xvel+i);
      // pass to op
      tensorop.setDomainBC({AMREX_D_DECL(mlmg_lobc[0],mlmg_lobc[1],mlmg_lobc[2])},
                           {AMREX_D_DECL(mlmg_hibc[0],mlmg_hibc[1],mlmg_hibc[2])});  
  
#if 0

  
      // set up level BCs
      {
        MultiFab crsedata;
        int ng = 0;

        if (level > 0) {
          auto& crse_ns = *(coarser->navier_stokes);
          crsedata.define(crse_ns.boxArray(), crse_ns.DistributionMap(), AMREX_SPACEDIM,
                          ng, MFInfo(),crse_ns.Factory());
          AmrLevel::FillPatch(crse_ns, crsedata, ng, cur_time, State_Type, Xvel,
                              AMREX_SPACEDIM);
          tensorop.setCoarseFineBC(&crsedata, crse_ratio[0]);
        }

        AmrLevel::FillPatch(*navier_stokes,Soln,soln_ng,cur_time,State_Type,Xvel,AMREX_SPACEDIM);

        tensorop.setLevelBC(0, &Soln);
      }

      {
        MultiFab acoef;
        acoef.define(grids, dmap, 1, 0, MFInfo(), navier_stokes->Factory());
        std::pair<Real,Real> scalars;
        // fixme? don't know why we're bothering with rhsscale....
        Real rhsscale = 1.0;
        const MultiFab& rho = (rho_flag == 1) ? rho_half : navier_stokes->rho_ctime;

        computeAlpha(acoef, scalars, a, b, rho,
                     1, &rhsscale, nullptr, 0,
                     nullptr, 0,
                     navier_stokes->Geom(),volume,parent->Geom(0).IsRZ());

        tensorop.setScalars(scalars.first, scalars.second);
        tensorop.setACoeffs(0, acoef);
      }



{
        Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
        for (int n = 0; n < BL_SPACEDIM; n++)
        {
          face_bcoef[n].define(area[n].boxArray(),area[n].DistributionMap(),1,0,MFInfo(),navier_stokes->Factory());
        }

        //      computeBeta(face_bcoef,betan,betaComp);
        computeBeta(face_bcoef,betan,betaComp,navier_stokes->Geom(),ap,
                    parent->Geom(0).IsRZ());
        tensorop.setShearViscosity(0, amrex::GetArrOfConstPtrs(face_bcoef));

#ifdef AMREX_USE_EB
        MultiFab cc_bcoef(grids,dmap,BL_SPACEDIM,0,MFInfo(),navier_stokes->Factory());
        EB_average_face_to_cellcenter(cc_bcoef, 0,
                                      amrex::GetArrOfConstPtrs(face_bcoef));
        tensorop.setEBShearViscosity(0, cc_bcoef);
#endif

        if (NavierStokesBase::S_in_vel_diffusion) {
          // remove the "divmusi" terms by setting kappa = (2/3) mu
          //
          Print()<<"WARNING: Hack to get rid of divU terms ...\n";
          Array<MultiFab,AMREX_SPACEDIM> kappa;
          Real twothirds = 2.0/3.0;
          for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
          {
            kappa[idim].define(face_bcoef[idim].boxArray(), face_bcoef[idim].DistributionMap(), 1, 0, MFInfo(),navier_stokes->Factory());
            MultiFab::Copy(kappa[idim], face_bcoef[idim], 0, 0, 1, 0);
            kappa[idim].mult(twothirds);
          }
          tensorop.setBulkViscosity(0, amrex::GetArrOfConstPtrs(kappa));
#ifdef AMREX_USE_EB
          cc_bcoef.mult(twothirds);
          tensorop.setEBBulkViscosity(0, cc_bcoef);
#endif
        }
      }  



 #endif 
  
  
  if (LES_model == "Smagorinsky") {
    amrex::Print() << "\n DEBUG WE DO SMAGORNISKY \n\n";
  }
  else
  {
    amrex::Abort("\n DEBUG DONT KNOW THIS LES MODEL \n\n");
  }
  
}
