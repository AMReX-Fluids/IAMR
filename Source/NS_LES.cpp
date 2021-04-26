

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
  NavierStokesBase& navier_stokes  = getLevel(level);
  const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes.Factory());
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

  //
  // Now that we have the gradients of velocity, we can compute the LES subgrid viscosity
  // 
      
  const auto dx = geom.CellSizeArray();

  for (MFIter mfi(Uvel,true); mfi.isValid(); ++mfi)
  {

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
  
      const Box& nbx = mfi.nodaltilebox(idim);
      Array4<Real      > dst = mu_LES[idim]->array(mfi);
      Array4<Real      > src = grad_Uvel[idim]->array(mfi);

      Real Cs_cst = LES_model == "Smagorinsky" ? smago_Cs_cst : sigma_Cs_cst;
      
      if (LES_model == "Smagorinsky") {
	
	if (ParallelDescriptor::IOProcessor() && getLESVerbose) {               
	  amrex::Print() << "\n in calc_mut_LES : WE DO SMAGORNISKY with constant " << Cs_cst << "\n\n";
	}

	AMREX_HOST_DEVICE_PARALLEL_FOR_4D (nbx, 1, i, j, k, n,
	{

            Real smag = 0;
            for (int i_symij = 0; i_symij < dim_fluxes; ++i_symij)
            {
              Real symij = src(i,j,k,i_symij) + src(i,j,k,i_symij);
              smag += symij * symij;
            }
                  
            smag = 0.5 * smag;
                  
            dst(i,j,k,n) = pow(Cs_cst * dx[idim],2) * sqrt(smag);
                    
	});
      } else if (LES_model == "Sigma") {

         //  Reference for the Sigma model
         //          Franck Nicoud, Hubert Baya Toda, Olivier Cabrit, Sanjeeb Bose, Jungil Lee
         //          Using singular values to build a subgrid-scale model for large eddy simulations
         //          Physics of Fluids, American Institute of Physics, 2011, 23 (8), pp.085106. ⟨10.1063/1.3623274⟩
         //          DOI : 10.1063/1.3623274


#if (AMREX_SPACEDIM < 3)
	amrex::Abort("FATAL ERROR in NS_LES.cpp: Sigma model is only for 3D");
#endif

	AMREX_HOST_DEVICE_PARALLEL_FOR_4D (nbx, 1, i, j, k, n,
        {

           Real G_11 = src(i,j,k,0)*src(i,j,k,0)  + src(i,j,k,1)*src(i,j,k,1)  +  src(i,j,k,2)*src(i,j,k,2);
           Real G_12 = src(i,j,k,0)*src(i,j,k,3)  + src(i,j,k,1)*src(i,j,k,4)  +  src(i,j,k,2)*src(i,j,k,5);
           Real G_13 = src(i,j,k,0)*src(i,j,k,6)  + src(i,j,k,1)*src(i,j,k,7)  +  src(i,j,k,2)*src(i,j,k,8);
           Real G_22 = src(i,j,k,3)*src(i,j,k,3)  + src(i,j,k,4)*src(i,j,k,4)  +  src(i,j,k,5)*src(i,j,k,5);
           Real G_23 = src(i,j,k,3)*src(i,j,k,6)  + src(i,j,k,4)*src(i,j,k,7)  +  src(i,j,k,5)*src(i,j,k,8);
           Real G_33 = src(i,j,k,6)*src(i,j,k,6)  + src(i,j,k,7)*src(i,j,k,7)  +  src(i,j,k,8)*src(i,j,k,8);

           //     First invariant (trace)
           Real I1 = G_11 + G_22 + G_33;
 
           //     Second invariant (0.5 * (tr(G)^2 - tr(G^2))
           Real I2 = G_11*G_22 - pow(G_12,2) + G_22*G_33 - pow(G_23,2) + G_11*G_33 - pow(G_13,2);

           //     Third invariant (determinant)
           Real I3 = G_11*(G_22*G_33 - G_23*G_23) - G_12*(G_33*G_12 - G_13*G_23) + G_13*(G_12*G_23 - G_13*G_22);

           //     Rotation angles
           Real alpha1 = max(0.,pow(I1/3,2) - I2/3);

           if ( alpha1 == 0. ){

              dst(i,j,k,n) = 0.;

           }
           else
           {

              Real alpha2 = pow(I1/3,3) -I1*I2/6 + I3/2;

              Real alphaArg = (alpha2*sqrt(1/alpha1))/alpha1;

              // Keeping AlphaArg between -1 and 1

              if ( alphaArg > 1. ){
                alphaArg = 1.;
              }
              else if ( alphaArg < -1. ){
                alphaArg = -1.;
              }
              Real alpha3 = acos(alphaArg)/3;

              //       Singular values
              
              // Ensuring that sigma1 >= sigma 2 >= sigma3 >=0 so that mu_LES is positive

              Real sigma1 = sqrt(max(0. ,I1/3 + 2 * sqrt(alpha1) * cos(alpha3)));
              Real sigma2 = sqrt(max(0. ,I1/3 - 2 * sqrt(alpha1) * cos(M_PI/3 + alpha3)));
              Real sigma3 = sqrt(max(0. ,I1/3 - 2 * sqrt(alpha1) * cos(M_PI/3 - alpha3)));
 
              Real verysmall = 1.e-24;
              sigma2 = max(sigma3,sigma2);
              sigma1 = max(sigma2,sigma1);
              sigma1 = max(verysmall,sigma1);


	      //       Compute the sigma operator
              dst(i,j,k,n) = pow(Cs_cst * dx[idim],2) * ((sigma3 * (sigma1-sigma2) * (sigma2-sigma3)) / pow(sigma1,2));

           }
	});

      } else {

	amrex::Abort("\n DEBUG DONT KNOW THIS LES MODEL \n\n");

      }
                    
    }
    
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
