
#include <NavierStokesBase.H>
#include <AMReX_BLFort.H>


using namespace amrex;

//
// Virtual access function for getting the forcing terms for the
// velocities and scalars.  The base version computes a buoyancy.
//
// NOTE: This function returns a rho weighted source term.
//
// For conservative (i.e. do_mom_diff=1, do_cons_trac=1), velocities
// are integrated according to the equation
//
//     ui_t + uj ui_j = S_ui        ===> tforces = rho S_ui
//
// and scalars psi where (psi = rho q = Scal) as
//
//     psi_t + (uj psi)_j = S_psi   ===> tforces = S_psi = rho S_q
//
// For non-conservative, this rho-weighted source term will get divided
// by rho in the predict_velocity, velocity_advection, scalar_advection,
// and advection_update routines.
//
// For temperature (which is always non-conservative), we evolve
//
//     dT/dt - U dot grad T = [del dot lambda grad T + S_T] / (rho*c_p)
//     ===> tforces =  S_T/c_p
//
//
// For user-defined forcing, this means
//   - For conservative variables, the force term computed here gets used
//     as-is
//   - For non-conservative variables, the force term computed here is
//     divided by rho before use
//

void
NavierStokesBase::getForce (FArrayBox&       force,
                            const Box&       bx,
                            int              scomp, // first component in force
                            int              ncomp, // number of components
                            const Real       time,
                            const FArrayBox& State, // state data, may not contain all components; e.g. may be velocities only
                            const FArrayBox& Aux,     // auxiliary data
                            int              auxScomp,// first component in Aux
                            const MFIter&    /*mfi*/)
{
   if (ParallelDescriptor::IOProcessor() && getForceVerbose)
   {
       const int*  f_lo     = force.loVect();
       const int*  f_hi     = force.hiVect();
       const int*  v_lo     = State.loVect();
       const int*  v_hi     = State.hiVect();
       const int*  s_lo     = Aux.loVect();
       const int*  s_hi     = Aux.hiVect();

       amrex::Print() << "NavierStokesBase::getForce(): Entered..." << std::endl
                      << "time      = " << time << std::endl
                      << "scomp     = " << scomp << std::endl
                      << "ncomp     = " << ncomp << std::endl
                      << "auxScomp = " << auxScomp << std::endl;

       if  (ncomp==1) amrex::Print() << "Doing only component " << scomp << std::endl;
       else if (scomp==0 && ncomp==AMREX_SPACEDIM) amrex::Print() << "Doing velocities only" << std::endl;
       else if (scomp>=AMREX_SPACEDIM) amrex::Print() << "Doing " << ncomp << " component(s) starting with component " << scomp << std::endl;

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
       amrex::Print() << "NavierStokesBase::getForce(): State Domain:" << std::endl;
       amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << ") - "
                      << "(" << v_hi[0] << "," << v_hi[1] << ")" << std::endl;
       amrex::Print() << "NavierStokesBase::getForce(): Aux Domain:" << std::endl;
       amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << ") - "
                      << "(" << s_hi[0] << "," << s_hi[1] << ")" << std::endl;
#endif

       // Compute min/max
       for (int n=0; n<ncomp; n++) {
           amrex::Print() << "State comp " << scomp+n << " min/max "
                          << State.min<RunOn::Gpu>(scomp+n) << " / "
                          << State.max<RunOn::Gpu>(scomp+n) << std::endl;
       }
       for (int n=auxScomp; n<Aux.nComp(); n++) {
           amrex::Print() << "aux comp " << n << " min/max "
                          << Aux.min<RunOn::Gpu>(n) << " / "
                          << Aux.max<RunOn::Gpu>(n) << std::endl;
       }
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
     auto const& frc = force.array(scomp);
     auto const& aux = Aux.array(auxScomp);
     const Real grav = gravity;

     if ( std::abs(grav) > 0.0001) {
       amrex::ParallelFor(bx, [frc, aux, grav]
       AMREX_GPU_DEVICE(int i, int j, int k) noexcept
       {
         frc(i,j,k,0) = Real(0.0);
#if ( AMREX_SPACEDIM == 2 )
         frc(i,j,k,1) = grav*aux(i,j,k,0);
#elif ( AMREX_SPACEDIM == 3 )
         frc(i,j,k,1) = Real(0.0);
         frc(i,j,k,2) = grav*aux(i,j,k,0);
#endif
       });
     }
     else {
       force.setVal<RunOn::Gpu>(0.0, bx, Xvel, AMREX_SPACEDIM);
     }
   }

   //
   // Scalar forcing
   //
   // During a regular timestep, getForce is called on scalars only.
   // During the multilevel sync, scalars are done with velocity.
   //
   int scomp_scal = -1;
   int ncomp_scal = -1;
   if ( scomp >= AMREX_SPACEDIM ) {
       // Doing only scalars
       scomp_scal = 0;
       ncomp_scal = ncomp;
   }
   // Recall that we will only get here if previous block is false,
   // i.e. if scomp < AMREX_SPACEDIM
   else if ( scomp+ncomp > AMREX_SPACEDIM) {
       // Doing scalars with vel
       scomp_scal = Density;
       ncomp_scal = ncomp-Density;
   }

   if (ncomp_scal > 0) {
       force.setVal<RunOn::Gpu>(0.0, bx, scomp_scal, ncomp_scal);
       //
       // Or create user-defined forcing.
       // Recall we compute a density-weighted forcing term.
       //
       // auto const& frc  = force.array(scomp_scal);
       // amrex::ParallelFor(bx, ncomp_scal, [frc]
       // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
       // {
       //          frc(i,j,k,n) = ;
       //          frc(i,j,k,n) *= rho;
       // });
   }

   if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
       // Compute min/max
       for (int n=0; n<ncomp; n++) {
           amrex::Print() << "Force comp " << scomp+n << " min/max "
                          << force.min<RunOn::Gpu>(scomp+n) << " / "
                          << force.max<RunOn::Gpu>(scomp+n) << std::endl;
       }

      amrex::Print() << "NavierStokesBase::getForce(): Leaving..."
                     << std::endl << "---" << std::endl;
   }
}
