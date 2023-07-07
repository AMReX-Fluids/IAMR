#ifndef IAMR_TurbulentForcing_params_H_
#define IAMR_TurbulentForcing_params_H_

namespace TurbulentForcing {

    // function to generate the parameters
    void init_turbulent_forcing (const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& problo, const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& probhi);

    extern AMREX_GPU_MANAGED int verbose;
    // factor by which to reduce sampling for faster performance
    extern AMREX_GPU_MANAGED int ff_factor;
    // make the forcing divergence free?
    extern AMREX_GPU_MANAGED bool div_free_force;
    // how many modes to use
    extern AMREX_GPU_MANAGED int nmodes;
    // don't use any modes below mode_start. We probably don't need this
    extern AMREX_GPU_MANAGED int mode_start;

    constexpr int array_size = 33;
    constexpr int num_fdarray = 17;
    // forcedata will contain num_fdarray arrays of size (0,0,0)(array_size-1,array_size-1,array_size-1)
    extern amrex::Real* forcedata;
}
#endif
