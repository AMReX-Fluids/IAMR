//fixme, for writesingle level plotfile
#include<AMReX_PlotFileUtil.H>
//

#include <NavierStokesBase.H>

#include <AMReX_VisMF.H>



using namespace amrex;

//--------------------------------------------------------------------
// "On the fly" time averaging of the velocity
//
//  Just put ns.avg_interval = something (>0) to activate the feature.
//  avg_interval controls the interval of time-steps when a solution is taken into account.
//  The average is dumped in each plotfile.
//
//  Do not forget to add "velocity_average" in amr.derive_plot_vars.
//
//  If "compute_fluctuations" is turned on, it is going to compute RMS of velocity fluctuations.
//  The good practice is to start computing RMS of fluctuations 
//  when the average of the velocity has reached convergence.
//---------------------------------------------------------------------

void
NavierStokesBase::time_average(amrex::Real&  time_avg, amrex::Real&  time_avg_fluct, amrex::Real&  dt_avg, const Real& dt_level)

{
  dt_avg = dt_avg + dt_level;

  if (parent->levelSteps(0)%avg_interval == 0)
  {
    MultiFab& Sstate = get_new_data(State_Type);
    MultiFab& Savg   = get_new_data(Average_Type); 
    MultiFab& Savg_old   = get_old_data(Average_Type);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Savg,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.tilebox();
       auto const& S_state = Sstate.array(mfi,Xvel);
       auto const& S_avg   = Savg.array(mfi);
       auto const& S_avg_old   = Savg_old.array(mfi);
       int loc_compute_fluctuations = compute_fluctuations; //NavierStokesBase class cannot be accessed directly fron device

       amrex::ParallelFor(bx, BL_SPACEDIM, [S_state, S_avg, S_avg_old, dt_avg, time_avg, time_avg_fluct, loc_compute_fluctuations]
       AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
          S_avg(i,j,k,n) = S_avg_old(i,j,k,n) + dt_avg * S_state(i,j,k,n);
          S_avg_old(i,j,k,n) = S_avg(i,j,k,n);

          if (loc_compute_fluctuations == 1){
            amrex::Real vel_prime = S_state(i,j,k,n) - (S_avg(i,j,k,n)/(time_avg + dt_avg));
            S_avg(i,j,k,n+BL_SPACEDIM) = S_avg_old(i,j,k,n+BL_SPACEDIM) + dt_avg * vel_prime * vel_prime;
          }
          else{
            S_avg(i,j,k,n+BL_SPACEDIM) = 0.;
          }
          S_avg_old(i,j,k,n+BL_SPACEDIM) = S_avg(i,j,k,n+BL_SPACEDIM);
       });
    }

    time_avg = time_avg + dt_avg;
    if (compute_fluctuations == 1){
      time_avg_fluct = time_avg_fluct + dt_avg;
    }else{
      time_avg_fluct = 0.;
    }

    dt_avg = 0;

  }
}

