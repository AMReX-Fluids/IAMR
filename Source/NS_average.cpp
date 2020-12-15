//fixme, for writesingle level plotfile
#include<AMReX_PlotFileUtil.H>
//

#include <NavierStokesBase.H>

#include <AMReX_VisMF.H>



using namespace amrex;



void
NavierStokesBase::time_average(bool flag_init, amrex::Real&  time_avg, amrex::Real&  dt_avg, const Real& dt_level)

{

  if (ParallelDescriptor::IOProcessor() ) {
    amrex::Print() << "\n  WE ARE IN THE NEW ROUTINE FOR TIME AVERAGE \n\n";
  }

  if (flag_init)
  {

    dt_avg   = 0;
    time_avg = 0;

    amrex::Print() << " " << std::endl;
    amrex::Print() << "DEBUG AVEGRAGE after INIT dt_avg = " << dt_avg << std::endl;
    amrex::Print() << "DEBUG AVEGRAGE after INIT time_avg = " << time_avg << std::endl;


MultiFab& Sstate = get_new_data(State_Type);
     MultiFab& Savg   = get_new_data(Average_Type);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif  
      for (MFIter mfi(Savg,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& S_state = Sstate.array(mfi,Xvel);
         auto const& S_avg   = Savg.array(mfi);

         amrex::ParallelFor(bx, BL_SPACEDIM, [S_state, S_avg, dt_avg]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            S_avg(i,j,k,n) = S_state(i,j,k,n);
         });
      }


  }
  else
  {

    amrex::Print() << " " << std::endl;
    amrex::Print() << "DEBUG AVEGRAGE dt_level = " << dt_level << std::endl;
    amrex::Print() << "DEBUG AVEGRAGE before dt_avg = " << dt_avg << std::endl;

    dt_avg = dt_avg + dt_level;
    amrex::Print() << "DEBUG AVEGRAGE after dt_avg = " << dt_avg << std::endl;

    int avg_interval = 1;
    if (parent->levelSteps(0)%avg_interval == 0)
    {
      amrex::Print() << "DEBUG HELLO FROM MODULO " << std::endl;
 
     MultiFab& Sstate = get_new_data(State_Type);
     MultiFab& Savg   = get_new_data(Average_Type); 

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Savg,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& S_state = Sstate.array(mfi,Xvel);
         auto const& S_avg   = Savg.array(mfi);

         amrex::ParallelFor(bx, BL_SPACEDIM, [S_state, S_avg, dt_avg]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            S_avg(i,j,k,n) = S_avg(i,j,k,n) + dt_avg * S_state(i,j,k,n);
         });
      }
   

  //      for (MFIter mfi(Savg,true); mfi.isValid(); ++mfi)
  //      {
  //       amrex::Print() << Savg[mfi];
  //      }
 
    }


  }
  

/*  
 
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
  // Now that we have the gradients of velocity, we can compute the LES subgrid viscosity
  // 
      
  const auto dx = geom.CellSizeArray();

  FArrayBox fab_tmp;
  for (MFIter mfi(Uvel,true); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();

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

  */
}

