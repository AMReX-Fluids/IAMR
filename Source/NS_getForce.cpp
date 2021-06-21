
#include <NavierStokesBase.H>
#include <AMReX_BLFort.H>


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
                            int              scomp,
                            int              ncomp,
                            const Real       time,
                            const FArrayBox& Vel,
                            const FArrayBox& Scal,
                            int              scalScomp,
                            const MFIter&    mfi)
{

   const Real* VelDataPtr  = Vel.dataPtr();
   const Real* ScalDataPtr = Scal.dataPtr(scalScomp);

   const Real  grav     = gravity;
   const int*  f_lo     = force.loVect();
   const int*  f_hi     = force.hiVect();
   const int*  v_lo     = Vel.loVect();
   const int*  v_hi     = Vel.hiVect();
   const int*  s_lo     = Scal.loVect();
   const int*  s_hi     = Scal.hiVect();

   if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
      amrex::Print() << "NavierStokesBase::getForce(): Entered..." << std::endl 
                     << "time      = " << time << std::endl
                     << "scomp     = " << scomp << std::endl
                     << "ncomp     = " << ncomp << std::endl
                     << "scalScomp = " << scalScomp << std::endl;

      if (scomp==0)
	if  (ncomp==3) amrex::Print() << "Doing velocities only" << std::endl;
	else           amrex::Print() << "Doing all components" << std::endl;
      else if (scomp==3)
	if  (ncomp==1) amrex::Print() << "Doing density only" << std::endl;
	else           amrex::Print() << "Doing all scalars" << std::endl;
      else if (scomp==4) amrex::Print() << "Doing tracer only" << std::endl;
      else               amrex::Print() << "Doing individual scalar" << std::endl;

      amrex::Print() << "NavierStokesBase::getForce(): Filling Force on box:"
		     << bx << std::endl;
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
   } //end if(getForceVerbose)

   //
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

     if ( std::abs(grav) > 0.0001) {
       amrex::ParallelFor(bx, [frc, scal, grav]
       AMREX_GPU_DEVICE(int i, int j, int k) noexcept
       {
	 frc(i,j,k,0) = Real(0.0);
#if ( AMREX_SPACEDIM == 2 )
         frc(i,j,k,1) = grav*scal(i,j,k,0);
#elif ( AMREX_SPACEDIM == 3 )
         frc(i,j,k,1) = Real(0.0);
         frc(i,j,k,2) = grav*scal(i,j,k,0);
#endif
       });
     }
     else {
       force.setVal<RunOn::Gpu>(0.0, bx, Xvel, AMREX_SPACEDIM);
       // amrex::ParallelFor(bx, AMREX_SPACEDIM, [frc]
       // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
       // {
       // 	 frc(i,j,k,n) = 0.0_rt;
       // });
     }
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
