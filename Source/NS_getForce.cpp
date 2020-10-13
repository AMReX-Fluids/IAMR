
#include <NavierStokesBase.H>
#include <AMReX_BLFort.H>
#include <PROB_NS_F.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

//
// Virtual access function for getting the forcing terms for the
// velocities and scalars.  The base version computes a buoyancy.
//
// As NavierStokesBase is currently implemented.  Velocities are integrated
// according to the equation
//
//     ui_t + uj ui_j = S_ui        ===> tforces = rho S_ui
//
// and scalars psi where (psi = rho q) as
//
//     psi_t + (uj psi)_j = S_psi   ===> tforces = S_psi = rho S_q
//
// q is a concentration.  This function returns a rho weighted
// source term, which requires a division by rho in the predict_velocity
// and velocity_advection routines.
//

void
NavierStokesBase::getForce (FArrayBox&       force,
                            const Box&       bx,
                            int              ngrow,
                            int              scomp,
                            int              ncomp,
                            const Real       time,
                            const FArrayBox& Vel,
                            const FArrayBox& Scal,
                            const FArrayBox& scalGrad,
                            int              scalScomp)
{

   const Real* VelDataPtr  = Vel.dataPtr();
   const Real* ScalDataPtr = Scal.dataPtr(scalScomp);
   const Real* ScalGradDataPtr = scalGrad.dataPtr();

   const Real* dx       = geom.CellSize();
   const Real  grav     = gravity;
   const int*  f_lo     = force.loVect();
   const int*  f_hi     = force.hiVect();
   const int*  v_lo     = Vel.loVect();
   const int*  v_hi     = Vel.hiVect();
   const int*  s_lo     = Scal.loVect();
   const int*  s_hi     = Scal.hiVect();
   const int*  sgrad_lo     = scalGrad.loVect();
   const int*  sgrad_hi     = scalGrad.hiVect();
   const int   nscal    = NUM_SCALARS;
      int tracInd = 0;
      // double sigma = 0.0;
      Vector<Real> sigma(nTrac); Vector<Real> sigmaS(nTrac);
      for(int n = 0; n < nTrac; n++)
      {
          sigma[n] = 0.0; sigmaS[n] = 0.0;
      }
      ParmParse pp("ns");
      if(ncomp == nscal)
      {
          tracInd = 1;
      }
      else
      {
          tracInd = Tracer;
      }

      //compute phasic-specific surface tension coefficients
      if(nTrac == 3)
      {
          pp.query("sigma12", sigma[0]);
          pp.query("sigma13", sigma[1]);
          pp.query("sigma23", sigma[2]);
          sigmaS[0] = (sigma[1]-sigma[2]+sigma[0])/2.0; //sigmaS1;
          sigmaS[1] = sigma[0] - sigmaS[0]; // sigmaS2;
          sigmaS[2] = sigma[1] - sigmaS[0]; //sigmaS3;
      }
      else if(nTrac == 2)
      {
          pp.query("sigma12", sigma[0]);
          sigmaS[0] = sigma[0];
          sigmaS[1] = 0.0;
      }
      else
      {
          pp.query("sigma12", sigma[0]);
          sigmaS[0] = sigma[0];

      }

   if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
      amrex::Print() << "NavierStokesBase::getForce(): Entered..." << std::endl
                     << "time      = " << time << std::endl
                     << "scomp     = " << scomp << std::endl
                     << "ncomp     = " << ncomp << std::endl
                     << "ngrow     = " << ngrow << std::endl
                     << "scalScomp = " << scalScomp << std::endl;

      if (scomp==0)
	if  (ncomp==3) amrex::Print() << "Doing velocities only" << std::endl;
	else           amrex::Print() << "Doing all components" << std::endl;
      else if (scomp==3)
	if  (ncomp==1) amrex::Print() << "Doing density only" << std::endl;
	else           amrex::Print() << "Doing all scalars" << std::endl;
      else if (scomp==4) amrex::Print() << "Doing tracer only" << std::endl;
      else               amrex::Print() << "Doing individual scalar" << std::endl;

#if (AMREX_SPACEDIM == 3)
      amrex::Print() << "NavierStokesBase::getForce(): Force Domain:" << std::endl;
      amrex::Print() << "(" << f_lo[0] << "," << f_lo[1] << "," << f_lo[2] << ") - "
                     << "(" << f_hi[0] << "," << f_hi[1] << "," << f_hi[2] << ")" << std::endl;
      amrex::Print() << "NavierStokesBase::getForce(): Vel Domain:" << std::endl;
      amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << "," << v_lo[2] << ") - "
                     << "(" << v_hi[0] << "," << v_hi[1] << "," << v_hi[2] << ")" << std::endl;
      amrex::Print() << "NavierStokesBase::getForce(): Scal Domain:" << std::endl;
      amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << "," << s_lo[2] << ") - "
                     << "(" << s_hi[0] << "," << s_hi[1] << "," << s_hi[2] << ")" << std::endl;
#else
      amrex::Print() << "NavierStokesBase::getForce(): Force Domain:" << std::endl;
      amrex::Print() << "(" << f_lo[0] << "," << f_lo[1] << ") - "
                     << "(" << f_hi[0] << "," << f_hi[1] << ")" << std::endl;
      amrex::Print() << "NavierStokesBase::getForce(): Vel Domain:" << std::endl;
      amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << ") - "
                     << "(" << v_hi[0] << "," << v_hi[1] << ")" << std::endl;
      amrex::Print() << "NavierStokesBase::getForce(): Scal Domain:" << std::endl;
      amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << ") - "
                     << "(" << s_hi[0] << "," << s_hi[1] << ")" << std::endl;
#endif

      Vector<Real> velmin(AMREX_SPACEDIM), velmax(AMREX_SPACEDIM);
      Vector<Real> scalmin(NUM_SCALARS), scalmax(NUM_SCALARS);
      Vector<Real> scalgradmin(3*nTrac), scalgradmax(3*nTrac);
      for (int n=0; n<AMREX_SPACEDIM; n++) {
          velmin[n]= 1.e234;
          velmax[n]=-1.e234;
      }
      int ix = v_hi[0]-v_lo[0]+1;
      int jx = v_hi[1]-v_lo[1]+1;
#if (AMREX_SPACEDIM == 3)
      int kx = v_hi[2]-v_lo[2]+1;
      for (int k=0; k<kx; k++) {
#endif
         for (int j=0; j<jx; j++) {
            for (int i=0; i<ix; i++) {
               for (int n=0; n<AMREX_SPACEDIM; n++) {
#if (AMREX_SPACEDIM == 3)
                  int cell = ((n*kx+k)*jx+j)*ix+i;
#else
                  int cell = (n*jx+j)*ix+i;
#endif
                  Real v = VelDataPtr[cell];
                  if (v<velmin[n]) velmin[n] = v;
                  if (v>velmax[n]) velmax[n] = v;
               }
            }
         }
#if (AMREX_SPACEDIM == 3)
      }
#endif
      for (int n=0; n<AMREX_SPACEDIM; n++)
         amrex::Print() << "Vel  " << n << " min/max "
                        << velmin[n] << " / " << velmax[n] << std::endl;

      for (int n=0; n<NUM_SCALARS; n++) {
         scalmin[n]= 1.e234;
         scalmax[n]=-1.e234;
      }
      ix = s_hi[0]-s_lo[0]+1;
      jx = s_hi[1]-s_lo[1]+1;
#if (AMREX_SPACEDIM == 3)
      kx = s_hi[2]-s_lo[2]+1;
      for (int k=0; k<kx; k++) {
#endif
         for (int j=0; j<jx; j++) {
            for (int i=0; i<ix; i++) {
               for (int n=0; n<NUM_SCALARS; n++) {
#if (AMREX_SPACEDIM == 3)
                  int cell = ((n*kx+k)*jx+j)*ix+i;
#else
                  int cell = (n*jx+j)*ix+i;
#endif
                  Real s = ScalDataPtr[cell];
                  if (s<scalmin[n]) scalmin[n] = s;
                  if (s>scalmax[n]) scalmax[n] = s;
               }
            }
         }
#if (AMREX_SPACEDIM == 3)
      }
#endif
  for (int n=0; n<NUM_SCALARS; n++)
     amrex::Print() << "Scal " << n << " min/max " << scalmin[n]
                    << " / " << scalmax[n] << std::endl;

                    for (int n=0; n<3*nTrac; n++) {
                       scalgradmin[n]= 1.e234;
                       scalgradmax[n]=-1.e234;
                    }
                    ix = sgrad_hi[0]-sgrad_lo[0]+1;
                    jx = sgrad_hi[1]-sgrad_lo[1]+1;
              #if (AMREX_SPACEDIM == 3)
                    kx = sgrad_hi[2]-sgrad_lo[2]+1;
                    for (int k=0; k<kx; k++) {
              #endif
                       for (int j=0; j<jx; j++) {
                          for (int i=0; i<ix; i++) {
                             for (int n=0; n<3*nTrac; n++) {
              #if (AMREX_SPACEDIM == 3)
                                int cell = ((n*kx+k)*jx+j)*ix+i;
              #else
                                int cell = (n*jx+j)*ix+i;
              #endif
                                Real sg = ScalGradDataPtr[cell];
                                if (sg<scalgradmin[n]) scalgradmin[n] = sg;
                                if (sg>scalgradmax[n]) scalgradmax[n] = sg;
                             }
                          }
                       }
              #if (AMREX_SPACEDIM == 3)
                    }
              #endif
                    for (int n=0; n<3*nTrac; n++)
                       amrex::Print() << "Scal grad " << n << " min/max " << scalgradmin[n]
                                      << " / " << scalgradmax[n] << std::endl;
   } //end if(getForceVerbose)

   RealBox gridloc = RealBox(bx,geom.CellSize(),geom.ProbLo());

   // Here's the meat
   //
   // Velocity forcing
   //
   if ( scomp<AMREX_SPACEDIM ){
     AMREX_ALWAYS_ASSERT(scomp==Xvel);
     AMREX_ALWAYS_ASSERT(ncomp>=AMREX_SPACEDIM);
   }

   if ( scomp==Xvel ){
     //
     // TODO: add some switch for user-supplied/problem-dependent forcing
     //
        auto const& frc  = force.array(scomp);
        auto const& scal = Scal.array(scalScomp);
        auto const& scalGradArr = scalGrad.array(0);

        // if ( std::abs(grav) > 0.0001) {
          amrex::ParallelFor(bx, [frc, scal, grav, scalGradArr, sigmaS]
          AMREX_GPU_DEVICE(int i, int j, int k) noexcept
          {
              double Fsvx; double Fsvy; double Fsvz;

              if(nTrac == 3)
              {

                  for(int n = 0; n < nTrac; n++)
                  {
                      Fsvx += sigmaS[n]*scalGradArr(i,j,k,3*n)*scalGradArr(i,j,k,3*n+2);
                      Fsvy += sigmaS[n]*scalGradArr(i,j,k,3*n+1)*scalGradArr(i,j,k,3*n+2);
                  }
              }
              else
              {
                  Fsvx = sigmaS[0]*scalGradArr(i,j,k,0)*scalGradArr(i,j,k,2);
                  Fsvy = sigmaS[0]*scalGradArr(i,j,k,1)*scalGradArr(i,j,k,2);
              }


        frc(i,j,k,0) = Fsvx;
     #if ( AMREX_SPACEDIM == 2 )
            frc(i,j,k,1) = grav*scal(i,j,k,0) + Fsvy;
     #elif ( AMREX_SPACEDIM == 3 )
            frc(i,j,k,1) = Fsvy;
            frc(i,j,k,2) = grav*scal(i,j,k,0);
     #endif
          });
        // }
        // else {
        //   // force.setVal<RunOn::Gpu>(0.0, bx, Xvel, AMREX_SPACEDIM);
        //   amrex::ParallelFor(bx, AMREX_SPACEDIM, [frc]
        //   AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
        //   {
        //   	 frc(i,j,k,n) = 0.0_rt;
        //   });
        // }
   }
   //
   // Scalar forcing
   //
   if ( scomp >= AMREX_SPACEDIM ) {
     // Doing only scalars
     force.setVal<RunOn::Gpu>(0.0, bx, 0, ncomp);
     // auto const& frc  = force.array();
     // amrex::ParallelFor(bx, ncomp, [frc]
     // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
     // {
     // 	 frc(i,j,k,n) = 0.0_rt;
     // });
   }
   else if ( scomp+ncomp > AMREX_SPACEDIM) {
     // Doing scalars with vel
     force.setVal<RunOn::Gpu>(0.0, bx, Density, ncomp-Density);
     // auto const& frc  = force.array(Density);
     // amrex::ParallelFor(bx, ncomp-Density, [frc]
     // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
     // {
     // 	 frc(i,j,k,n) = 0.0_rt;
     // });
   }

   if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
      Vector<Real> forcemin(ncomp);
      Vector<Real> forcemax(ncomp);
      for (int n=0; n<ncomp; n++) {
         forcemin[n]= 1.e234;
         forcemax[n]=-1.e234;
      }

      int ix = f_hi[0]-f_lo[0]+1;
      int jx = f_hi[1]-f_lo[1]+1;
#if (AMREX_SPACEDIM == 3)
      int kx = f_hi[2]-f_lo[2]+1;
      for (int k=0; k<kx; k++) {
#endif
         for (int j=0; j<jx; j++) {
            for (int i=0; i<ix; i++) {
               for (int n=0; n<ncomp; n++) {
#if (AMREX_SPACEDIM == 3)
                  int cell = ((n*kx+k)*jx+j)*ix+i;
#else
                  int cell = (n*jx+j)*ix+i;
#endif
                  Real f = force.dataPtr()[cell];
                  if (f<forcemin[n]) forcemin[n] = f;
                  if (f>forcemax[n]) forcemax[n] = f;
               }
            }
         }
#if (AMREX_SPACEDIM == 3)
      }
#endif
      for (int n=0; n<ncomp; n++)
         amrex::Print() << "Force " << n+scomp << " min/max " << forcemin[n]
                        << " / " << forcemax[n] << std::endl;

      amrex::Print() << "NavierStokesBase::getForce(): Leaving..."
                     << std::endl << "---" << std::endl;
   }
}


// void
// NavierStokesBase::getForce (FArrayBox&       force,
//                             const Box&       bx,
//                             int              ngrow,
//                             int              scomp,
//                             int              ncomp,
//                             const Real       time,
//                             const FArrayBox& Vel,
//                             const FArrayBox& Scal,
//                             const FArrayBox& scalGrad,
//                             int              scalScomp)
// {
//
//    const Real* VelDataPtr  = Vel.dataPtr();
//    const Real* ScalDataPtr = Scal.dataPtr(scalScomp);
//    const Real* ScalGradDataPtr = scalGrad.dataPtr(scalScomp);
//
//    const Real* dx       = geom.CellSize();
//    const Real  grav     = gravity;
//    const int*  f_lo     = force.loVect();
//    const int*  f_hi     = force.hiVect();
//    const int*  v_lo     = Vel.loVect();
//    const int*  v_hi     = Vel.hiVect();
//    const int*  s_lo     = Scal.loVect();
//    const int*  s_hi     = Scal.hiVect();
//    const int*  sgrad_lo     = scalGrad.loVect();
//    const int*  sgrad_hi     = scalGrad.hiVect();
//    const int   nscal    = NUM_SCALARS;
//    int tracInd = 0;
//    // double sigma = 0.0;
//    Vector<Real> sigma(nTrac); Vector<Real> sigmaS(nTrac);
//    for(int n = 0; n < nTrac; n++)
//    {
//        sigma[n] = 0.0; sigmaS[n] = 0.0;
//    }
//    ParmParse pp("ns");
//    if(ncomp == nscal)
//    {
//        tracInd = 1;
//    }
//    else
//    {
//        tracInd = Tracer;
//    }
//
//    //compute phasic-specific surface tension coefficients
//    if(nTrac == 3)
//    {
//        pp.query("sigma12", sigma[0]);
//        pp.query("sigma13", sigma[1]);
//        pp.query("sigma23", sigma[2]);
//        sigmaS[0] = (sigma[1]-sigma[2]+sigma[0])/2.0; //sigmaS1;
//        sigmaS[1] = sigma[0] - sigmaS[0]; // sigmaS2;
//        sigmaS[2] = sigma[1] - sigmaS[0]; //sigmaS3;
//    }
//    else if(nTrac == 2)
//    {
//        pp.query("sigma12", sigma[0]);
//        sigmaS[0] = sigma[0];
//        sigmaS[1] = 0.0;
//    }
//    else
//    {
//        pp.query("sigma12", sigma[0]);
//        sigmaS[0] = sigma[0];
//
//    }
//
//    // amrex::Print() << "NavierStokesBase::getForce(): sigma:" << sigmaS[0] <<  std::endl;
//    if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
//       amrex::Print() << "NavierStokesBase::getForce(): Entered..." << std::endl
//                      << "time      = " << time << std::endl
//                      << "scomp     = " << scomp << std::endl
//                      << "ncomp     = " << ncomp << std::endl
//                      << "ngrow     = " << ngrow << std::endl
//                      << "scalScomp = " << scalScomp << std::endl;
//
//       if (scomp==0)
// 	if  (ncomp==3) amrex::Print() << "Doing velocities only" << std::endl;
// 	else           amrex::Print() << "Doing all components" << std::endl;
//       else if (scomp==3)
// 	if  (ncomp==1) amrex::Print() << "Doing density only" << std::endl;
// 	else           amrex::Print() << "Doing all scalars" << std::endl;
//       else if (scomp==4) amrex::Print() << "Doing tracer only" << std::endl;
//       else               amrex::Print() << "Doing individual scalar" << std::endl;
//
// #if (AMREX_SPACEDIM == 3)
//       amrex::Print() << "NavierStokesBase::getForce(): Force Domain:" << std::endl;
//       amrex::Print() << "(" << f_lo[0] << "," << f_lo[1] << "," << f_lo[2] << ") - "
//                      << "(" << f_hi[0] << "," << f_hi[1] << "," << f_hi[2] << ")" << std::endl;
//       amrex::Print() << "NavierStokesBase::getForce(): Vel Domain:" << std::endl;
//       amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << "," << v_lo[2] << ") - "
//                      << "(" << v_hi[0] << "," << v_hi[1] << "," << v_hi[2] << ")" << std::endl;
//       amrex::Print() << "NavierStokesBase::getForce(): Scal Domain:" << std::endl;
//       amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << "," << s_lo[2] << ") - "
//                      << "(" << s_hi[0] << "," << s_hi[1] << "," << s_hi[2] << ")" << std::endl;
// #else
//       amrex::Print() << "NavierStokesBase::getForce(): Force Domain:" << std::endl;
//       amrex::Print() << "(" << f_lo[0] << "," << f_lo[1] << ") - "
//                      << "(" << f_hi[0] << "," << f_hi[1] << ")" << std::endl;
//       amrex::Print() << "NavierStokesBase::getForce(): Vel Domain:" << std::endl;
//       amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << ") - "
//                      << "(" << v_hi[0] << "," << v_hi[1] << ")" << std::endl;
//       amrex::Print() << "NavierStokesBase::getForce(): Scal Domain:" << std::endl;
//       amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << ") - "
//                      << "(" << s_hi[0] << "," << s_hi[1] << ")" << std::endl;
// #endif
//
//       Vector<Real> velmin(AMREX_SPACEDIM), velmax(AMREX_SPACEDIM);
//       Vector<Real> scalmin(NUM_SCALARS), scalmax(NUM_SCALARS);
//       for (int n=0; n<AMREX_SPACEDIM; n++) {
//           velmin[n]= 1.e234;
//           velmax[n]=-1.e234;
//       }
//       int ix = v_hi[0]-v_lo[0]+1;
//       int jx = v_hi[1]-v_lo[1]+1;
// #if (AMREX_SPACEDIM == 3)
//       int kx = v_hi[2]-v_lo[2]+1;
//       for (int k=0; k<kx; k++) {
// #endif
//          for (int j=0; j<jx; j++) {
//             for (int i=0; i<ix; i++) {
//                for (int n=0; n<AMREX_SPACEDIM; n++) {
// #if (AMREX_SPACEDIM == 3)
//                   int cell = ((n*kx+k)*jx+j)*ix+i;
// #else
//                   int cell = (n*jx+j)*ix+i;
// #endif
//                   Real v = VelDataPtr[cell];
//                   if (v<velmin[n]) velmin[n] = v;
//                   if (v>velmax[n]) velmax[n] = v;
//                }
//             }
//          }
// #if (AMREX_SPACEDIM == 3)
//       }
// #endif
//       for (int n=0; n<AMREX_SPACEDIM; n++)
//          amrex::Print() << "Vel  " << n << " min/max "
//                         << velmin[n] << " / " << velmax[n] << std::endl;
//
//       for (int n=0; n<NUM_SCALARS; n++) {
//          scalmin[n]= 1.e234;
//          scalmax[n]=-1.e234;
//       }
//       ix = s_hi[0]-s_lo[0]+1;
//       jx = s_hi[1]-s_lo[1]+1;
// #if (AMREX_SPACEDIM == 3)
//       kx = s_hi[2]-s_lo[2]+1;
//       for (int k=0; k<kx; k++) {
// #endif
//          for (int j=0; j<jx; j++) {
//             for (int i=0; i<ix; i++) {
//                for (int n=0; n<NUM_SCALARS; n++) {
// #if (AMREX_SPACEDIM == 3)
//                   int cell = ((n*kx+k)*jx+j)*ix+i;
// #else
//                   int cell = (n*jx+j)*ix+i;
// #endif
//                   Real s = ScalDataPtr[cell];
//                   if (s<scalmin[n]) scalmin[n] = s;
//                   if (s>scalmax[n]) scalmax[n] = s;
//                }
//             }
//          }
// #if (AMREX_SPACEDIM == 3)
//       }
// #endif
//       for (int n=0; n<NUM_SCALARS; n++)
//          amrex::Print() << "Scal " << n << " min/max " << scalmin[n]
//                         << " / " << scalmax[n] << std::endl;
//    } //end if(getForceVerbose)
//
//
//    RealBox gridloc = RealBox(bx,geom.CellSize(),geom.ProbLo());
//    const Real* xlo = gridloc.lo();
//    const Real* xhi = gridloc.hi();
//    // Here's the meat
//    // FORT_MAKEFORCE (&time,
//    //                 BL_TO_FORTRAN_ANYD(force),
//    //                 BL_TO_FORTRAN_ANYD(Vel),
//    //                 BL_TO_FORTRAN_N_ANYD(Scal,scalScomp),
//    //                 // BL_TO_FORTRAN_ANYD(scalGrad),
//    //                 dx,
//    //                 gridloc.lo(),
//    //                 gridloc.hi(),
//    //                 &grav, &sigma,&scomp,&ncomp,&nscal,&getForceVerbose);
//
//  // Array4<Real> const& forceArray = force.array();
//  // //std::cout <<  forceArray.nComp() << std::endl;
//  // Array4<const Real> scalArray = Scal.array();
//  // Array4<const Real> scalGradArr = scalGrad.array();
//  double rho = 0.0;
//  double trac = 0.0; double trac2 = 0.0; double trac3 = 0.0;
//  Real dTdx, dTdy, d2Tdx2, d2Tdy2, d2Tdydx, magn, curv;
//  Vector<Real> scalgradmax(3); Vector<Real> scalgradmin(3);
//  Vector<Real> tempmax(3); Vector<Real> tempmin(3);
//  for(int i = 0; i < 3; i++)
//  {
//      tempmax[i] = 0.0; tempmin[i] = 0.0; scalgradmax[i] = 0.0; scalgradmin[i] = 0.0;
//  }
//
//  if ( scomp==Xvel ){
//    //
//    // TODO: add some switch for user-supplied/problem-dependent forcing
//    //
//    auto const& frc  = force.array(scomp);
//    auto const& scal = Scal.array(scalScomp);
//    auto const& scalGradArr = scalGrad.array();
//
//    if ( std::abs(grav) > 0.0001) {
//      amrex::ParallelFor(bx, [frc, scal, grav, scalGradArr, sigmaS]
//      AMREX_GPU_DEVICE(int i, int j, int k) noexcept
//      {
//          double Fsvx = 0.0; double Fsvy = 0.0; double Fsvz = 0.0;
//
//          if(nTrac == 3)
//          {
//
//              for(int n = 0; n < nTrac; n++)
//              {
//                  Fsvx += sigmaS[n]*scalGradArr(i,j,k,3*n)*scalGradArr(i,j,k,3*n+2);
//                  Fsvy += sigmaS[n]*scalGradArr(i,j,k,3*n+1)*scalGradArr(i,j,k,3*n+2);
//              }
//          }
//          else
//          {
//              Fsvx = sigmaS[0]*scalGradArr(i,j,k,0)*scalGradArr(i,j,k,2);
//              Fsvy = sigmaS[0]*scalGradArr(i,j,k,1)*scalGradArr(i,j,k,2);
//          }
//
//
//    frc(i,j,k,0) = Fsvx;
// #if ( AMREX_SPACEDIM == 2 )
//        frc(i,j,k,1) = grav*scal(i,j,k,0) + Fsvy;
// #elif ( AMREX_SPACEDIM == 3 )
//        frc(i,j,k,1) = Fsvy;
//        frc(i,j,k,2) = grav*scal(i,j,k,0);
// #endif
//      });
//    }
//    else {
//      force.setVal<RunOn::Gpu>(0.0, bx, Xvel, AMREX_SPACEDIM);
//      // amrex::ParallelFor(bx, AMREX_SPACEDIM, [frc]
//      // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
//      // {
//      // 	 frc(i,j,k,n) = 0.0_rt;
//      // });
//    }
//  }
//  //
//  // Scalar forcing
//  //
//  if ( scomp >= AMREX_SPACEDIM ) {
//    // Doing only scalars
//    force.setVal<RunOn::Gpu>(0.0, bx, 0, ncomp);
//    // auto const& frc  = force.array();
//    // amrex::ParallelFor(bx, ncomp, [frc]
//    // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
//    // {
//    // 	 frc(i,j,k,n) = 0.0_rt;
//    // });
//  }
//  else if ( scomp+ncomp > AMREX_SPACEDIM) {
//    // Doing scalars with vel
//    force.setVal<RunOn::Gpu>(0.0, bx, Density, ncomp-Density);
//    // auto const& frc  = force.array(Density);
//    // amrex::ParallelFor(bx, ncomp-Density, [frc]
//    // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
//    // {
//    // 	 frc(i,j,k,n) = 0.0_rt;
//    // });
//  }
//
//  if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
//     Vector<Real> forcemin(ncomp);
//     Vector<Real> forcemax(ncomp);
//     for (int n=0; n<ncomp; n++) {
//        forcemin[n]= 1.e234;
//        forcemax[n]=-1.e234;
//     }
//
//     int ix = f_hi[0]-f_lo[0]+1;
//     int jx = f_hi[1]-f_lo[1]+1;
// #if (AMREX_SPACEDIM == 3)
//     int kx = f_hi[2]-f_lo[2]+1;
//     for (int k=0; k<kx; k++) {
// #endif
//        for (int j=0; j<jx; j++) {
//           for (int i=0; i<ix; i++) {
//              for (int n=0; n<ncomp; n++) {
// #if (AMREX_SPACEDIM == 3)
//                 int cell = ((n*kx+k)*jx+j)*ix+i;
// #else
//                 int cell = (n*jx+j)*ix+i;
// #endif
//                 Real f = force.dataPtr()[cell];
//                 if (f<forcemin[n]) forcemin[n] = f;
//                 if (f>forcemax[n]) forcemax[n] = f;
//              }
//           }
//        }
// #if (AMREX_SPACEDIM == 3)
//     }
// #endif
//     for (int n=0; n<ncomp; n++)
//        amrex::Print() << "Force " << n+scomp << " min/max " << forcemin[n]
//                       << " / " << forcemax[n] << std::endl;
//
//     amrex::Print() << "NavierStokesBase::getForce(): Leaving..."
//                    << std::endl << "---" << std::endl;
//  }
//
//
//
//
// //  if(scomp == 0)
// //  {
// //      for(int i = f_lo[0]; i <=f_hi[0]; i++)
// //      {
// //          for(int j = f_lo[1]; j <=f_hi[1]; j++)
// //          {
// //              #if (AMREX_SPACEDIM == 3)
// //             for(int k = f_lo[2]; k <= f_hi[2]; k++)
// //             {
// //             #else
// //             int k = 0;
// //             #endif
// //
// //             if(scalScomp == 0)
// //             {
// //                 rho = scalArray(i,j,k,0);
// //                 trac = scalArray(i,j,k,1);
// //             }
// //             else
// //             {
// //                 rho = scalArray(i,j,k,Density);
// //                 trac = scalArray(i,j,k,Tracer);
// //             }
// //
// //
// //             // if(isnan(scalGradArr(i,j,0,1)) || isnan(scalGradArr(i,j,0,0)))
// //             // {
// //             //         std::cout << "Scal grad = nan!" << std::endl;
// //             //         std::cout << i << " " << j << std::endl;
// //             //         exit(1);
// //             // }
// //              // if(dTdx == 0 && dTdy == 0)
// //              // {
// //              for(int n = 0; n < 3; n++)
// //              {
// //                  if(tempmax[n] < scalGradArr(i,j,k,n))
// //                  {
// //                      tempmax[n] = scalGradArr(i,j,k,n);
// //                  }
// //                  if(tempmin[n] > scalGradArr(i,j,k,n))
// //                  {
// //                      tempmin[n] = scalGradArr(i,j,k,n);
// //                  }
// //              }
// //
// //              // }
// //                 if(nTrac == 3)
// //                 {
// //
// //                     for(int n = 0; n < nTrac; n++)
// //                     {
// //                         Fsvx += sigmaS[n]*scalGradArr(i,j,k,3*n)*scalGradArr(i,j,k,3*n+2);
// //                         Fsvy += sigmaS[n]*scalGradArr(i,j,k,3*n+1)*scalGradArr(i,j,k,3*n+2);
// //                     }
// //                 }
// //                 else
// //                 {
// //                     Fsvx = sigmaS[0]*scalGradArr(i,j,k,0)*scalGradArr(i,j,k,2);
// //                     Fsvy = sigmaS[0]*scalGradArr(i,j,k,1)*scalGradArr(i,j,k,2);
// //                 }
// //
// //                  // Fsvx = 2.0*sigmaS[0]*trac*scalGradArr(i,j,0,0)*scalGradArr(i,j,0,2);
// //                  // Fsvy = 2.0*sigmaS[0]*trac*scalGradArr(i,j,0,1)*scalGradArr(i,j,0,2);
// //
// //                  // Fsvx = 2.0*sigma*trac*dTdx*curv;
// //                  // Fsvy = 2.0*sigma*trac*dTdy*curv;
// //
// //                  forceArray(i,j,k,Xvel) =  Fsvx;
// //                  forceArray(i,j,k,Yvel) = rho*grav + Fsvy;
// //                  #if (AMREX_SPACEDIM == 3)
// //                  forceArray(i,j,k,Zvel) =  0.0;
// //                  #endif
// //                   // if(i == 16 && j == 16)
// //                   // {
// //                   //     std::cout << "dTdx = " << scalGradArr(i,j,0, 0) << std::endl;
// //                   //
// //                   // }
// //                   #if (AMREX_SPACEDIM == 3)
// //                     }
// //                 #endif
// //
// //          }
// //      }
// //      if(ParallelDescriptor::IOProcessor() && getForceVerbose)
// //      {
// //
// //
// //      for(int n = 0; n < 3; n++)
// //      {
// //          scalgradmax[n] = tempmax[n];
// //          scalgradmin[n] = tempmin[n];
// //          std::cout << "scal grad max " << n << " = " << scalgradmax[n] << std::endl;
// //          std::cout << "scal grad min " << n << " = " << scalgradmin[n] << std::endl;
// //      }
// //      }
// //  }
// //
// //
// //  if((scomp+ncomp) > AMREX_SPACEDIM)
// //  {
// //      // for(int n = 0; n < NUM_SCALARS; n++)
// //      // {
// //          for(int i = f_lo[0]; i <= f_hi[0]; i++)
// //          {
// //              for(int j = f_lo[1]; j <= f_hi[1]; j++)
// //              {
// //                  #if (AMREX_SPACEDIM == 3)
// //
// //                  for(int k = f_lo[2]; k <= f_hi[2]; k++)
// //                  {
// //                 #else
// //                     int k = 0;
// //                 #endif
// //                  // if(n == AMREX_SPACEDIM)
// //                  // {
// //                     //std::cout << AMREX_SPACEDIM << " " << scomp << " " << ncomp << "   ";
// //                     if(ncomp == 1)
// //                     {
// //                         forceArray(i,j,k,0) = 0.0;
// //                     }
// //                     else
// //                     {
// //                         for(int n = 0; n < ncomp; n++)
// //                         {
// //                             forceArray(i,j,k,n) = 0.0;
// //                         }
// //
// //                         // std::cout << forceArray(i,j,0,1) << " ";
// //                     }
// //                  // }
// //                  // else if(n == AMREX_SPACEDIM + 1)
// //                  // {
// //                  //     forceArray(i,j,0,n) = 0.0;
// //                  // }
// //                  // else
// //                  // {
// //                  //     forceArray(i,j,0,n) = 0.0;
// //                  // }
// //              }
// //              #if (AMREX_SPACEDIM == 3)
// //              }
// //              #endif
// //          }
// //
// //      // }
// //  }
// //
// // //      if ( std::abs(grav) > 0.0001) {
// // //        amrex::ParallelFor(bx, [frc, scal, grav]
// // //        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
// // //        {
// // // 	 frc(i,j,k,0) = 0.0_rt;
// // // #if ( AMREX_SPACEDIM == 2 )
// // //          frc(i,j,k,1) = grav*scal(i,j,k,0);
// // // #elif ( AMREX_SPACEDIM == 3 )
// // //          frc(i,j,k,1) = 0.0_rt;
// // //          frc(i,j,k,2) = grav*scal(i,j,k,0);
// // // #endif
// // //        });
// // //      }
// // //      else {
// // //        force.setVal<RunOn::Gpu>(0.0, bx, Xvel, AMREX_SPACEDIM);
// // //        // amrex::ParallelFor(bx, AMREX_SPACEDIM, [frc]
// // //        // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
// // //        // {
// // //        // 	 frc(i,j,k,n) = 0.0_rt;
// // //        // });
// // //      }
// // //    }
// // //    //
// // //    // Scalar forcing
// // //    //
// // //    if ( scomp >= AMREX_SPACEDIM ) {
// // //      // Doing only scalars
// // //      force.setVal<RunOn::Gpu>(0.0, bx, 0, ncomp);
// // //      // auto const& frc  = force.array();
// // //      // amrex::ParallelFor(bx, ncomp, [frc]
// // //      // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
// // //      // {
// // //      // 	 frc(i,j,k,n) = 0.0_rt;
// // //      // });
// // //    }
// // //    else if ( scomp+ncomp > AMREX_SPACEDIM) {
// // //      // Doing scalars with vel
// // //      force.setVal<RunOn::Gpu>(0.0, bx, Density, ncomp-Density);
// // //      // auto const& frc  = force.array(Density);
// // //      // amrex::ParallelFor(bx, ncomp-Density, [frc]
// // //      // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
// // //      // {
// // //      // 	 frc(i,j,k,n) = 0.0_rt;
// // //      // });
// // //    }
// //
// //    if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
// //       Vector<Real> forcemin(ncomp);
// //       Vector<Real> forcemax(ncomp);
// //       for (int n=0; n<ncomp; n++) {
// //          forcemin[n]= 1.e234;
// //          forcemax[n]=-1.e234;
// //       }
// //
// //       int ix = f_hi[0]-f_lo[0]+1;
// //       int jx = f_hi[1]-f_lo[1]+1;
// // #if (AMREX_SPACEDIM == 3)
// //       int kx = f_hi[2]-f_lo[2]+1;
// //       for (int k=0; k<kx; k++) {
// // #endif
// //          for (int j=0; j<jx; j++) {
// //             for (int i=0; i<ix; i++) {
// //                for (int n=0; n<ncomp; n++) {
// // #if (AMREX_SPACEDIM == 3)
// //                   int cell = ((n*kx+k)*jx+j)*ix+i;
// // #else
// //                   int cell = (n*jx+j)*ix+i;
// // #endif
// //                   Real f = force.dataPtr()[cell];
// //                   if (f<forcemin[n]) forcemin[n] = f;
// //                   if (f>forcemax[n]) forcemax[n] = f;
// //                }
// //             }
// //          }
// // #if (AMREX_SPACEDIM == 3)
// //       }
// // #endif
// //       for (int n=0; n<ncomp; n++)
// //          amrex::Print() << "Force " << n+scomp << " min/max " << forcemin[n]
// //                         << " / " << forcemax[n] << std::endl;
// //
// //       amrex::Print() << "NavierStokesBase::getForce(): Leaving..."
// //                      << std::endl << "---" << std::endl;
// //    }
// // }
