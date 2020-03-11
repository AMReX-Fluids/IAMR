

//fixme, for writesingle level plotfile
#include<AMReX_PlotFileUtil.H>
//

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
NavierStokesBase::calc_mut_LES(MultiFab* mu_LES[BL_SPACEDIM], const Real time)

{

  if (ParallelDescriptor::IOProcessor() && getLESVerbose) {
    amrex::Print() << "\n in calc_mut_LES : WE ARE IN THE NEW ROUTINE FOR LES \n\n";
  }
  
  //
  // We get the state data at the current time in order to get the velocity
  //
  
  auto whichTime = which_time(State_Type,time);
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  MultiFab& Sstate = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);
  
  int nGrow = 1;
  FillPatchIterator fpi(*this,Sstate,nGrow,time,State_Type,Xvel,AMREX_SPACEDIM);
  MultiFab& Uvel=fpi.get_mf();
  
  //
  // Creating the TensorOp object to compute gradients of velocity at each face
  //
  
  LPInfo info;

#ifdef AMREX_USE_EB
  const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory());
  MLEBTensorOp tensorop({geom}, {grids}, {dmap}, info, {ebf});

  std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM>areafrac = ebf->getAreaFrac();
#else
  MLTensorOp tensorop({geom}, {grids}, {dmap}, info);
#endif

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
  {
    MultiFab crsedata;
    int ng = 0;
    const int soln_ng = 1;
      
    if (level > 0) {
      //auto& crse_ns = *(coarser->navier_stokes);
      NavierStokesBase& crse_ns  = getLevel(level-1);
      crsedata.define(crse_ns.boxArray(), crse_ns.DistributionMap(), AMREX_SPACEDIM,
                      ng, MFInfo(),crse_ns.Factory());
      AmrLevel::FillPatch(crse_ns, crsedata, ng, time, State_Type, Xvel,
                          AMREX_SPACEDIM);
      tensorop.setCoarseFineBC(&crsedata, crse_ratio[0]);
     }
      
     AmrLevel::FillPatch(*this,Uvel,soln_ng,time,State_Type,Xvel,AMREX_SPACEDIM);
      
     tensorop.setLevelBC(0, &Uvel);
  }
              
  const int dim_fluxes = pow(AMREX_SPACEDIM,2);
  
  FluxBoxes fb(this,dim_fluxes);
  MultiFab** tensorflux = fb.get();
  std::array<MultiFab*,AMREX_SPACEDIM> grad_Uvel{AMREX_D_DECL(tensorflux[0], tensorflux[1], tensorflux[2])};
 
  tensorop.compVelGrad(0,{grad_Uvel},{Uvel},MLLinOp::Location::FaceCenter);

//VisMF::Write(*grad_Uvel[0],"fluxes_x");
//VisMF::Write(*grad_Uvel[1],"fluxes_y");
  
//static int count_LES_src=0;
//count_LES_src++;
//amrex::WriteSingleLevelPlotfile("grad_Uvel_"+std::to_string(count_LES_src), *grad_Uvel[0], {"comp0","comp1","comp2","comp3"}, geom, 0.0, 0);

  
 
 
  //
  // Now that we have the gradients of velocity, we can compute the LES subgrid viscosity
  // 
      
  const auto dx = geom.CellSizeArray();

// amrex::Print() << "\n DEBUG geom.CellSizeArray() " << dx[0]  << "\n\n";

  FArrayBox fab_tmp;
  for (MFIter mfi(Uvel,true); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
  
  //amrex::Print() << "\n WE ARE IN IDIM " << idim << "\n ";
  
      const Box& nbx = mfi.nodaltilebox(idim);
      Array4<Real      > dst = mu_LES[idim]->array(mfi);
      Array4<Real      > src = grad_Uvel[idim]->array(mfi);

      AMREX_HOST_DEVICE_PARALLEL_FOR_4D (nbx, 1, i, j, k, n,
      {
         if (LES_model == "Smagorinsky") {
          
            if (ParallelDescriptor::IOProcessor() && getLESVerbose) {               
                amrex::Print() << "\n in calc_mut_LES : WE DO SMAGORNISKY with constant " << smago_Cs_cst << "\n\n";
            }
          
            Real smag = 0;
            for (int i_symij = 0; i_symij < dim_fluxes; ++i_symij)
            {
              Real symij = src(i,j,k,i_symij) + src(i,j,k,i_symij);
              smag += symij * symij;
            }
                  
            smag = 0.5 * smag;
                  
            dst(i,j,k,n) = pow(smago_Cs_cst * dx[idim],2) * sqrt(smag);
                    
         }
         else
         {
           amrex::Abort("\n DEBUG DONT KNOW THIS LES MODEL \n\n");
         }
                    
         //amrex::Print() << "\n MU_LES " << i << " " << j << " " << k << " " << dst(i,j,k,n);
      });
    }
    
  }

//static int count_LES=0;
//count_LES++;
//amrex::WriteSingleLevelPlotfile("mu_LES_"+std::to_string(count_LES), *mu_LES[0], {"mu_LES_x"}, geom, 0.0, 0);
//


//VisMF::Write(*mu_LES[0],"mu_LES_x");
//VisMF::Write(*mu_LES[1],"mu_LES_y");

 
  
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
