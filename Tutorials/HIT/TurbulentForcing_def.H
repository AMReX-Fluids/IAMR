#include <TurbulentForcing_params.H>
#include <depRand.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Arena.H>
#include <iamr_constants.H>

// factor by which to reduce sampling for faster performance
AMREX_GPU_MANAGED int TurbulentForcing::ff_factor;
// make the forcing divergence free?
AMREX_GPU_MANAGED bool TurbulentForcing::div_free_force;
// how many modes to use
AMREX_GPU_MANAGED int TurbulentForcing::nmodes;
// don't use any modes below mode_start. We probably don't need this
AMREX_GPU_MANAGED int TurbulentForcing::mode_start;
// Diagnostic print outs
AMREX_GPU_MANAGED int TurbulentForcing::verbose;

amrex::Real* TurbulentForcing::forcedata;


void
TurbulentForcing::init_turbulent_forcing (const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& problo, const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& probhi)
{
    using namespace amrex;

    // Start with checks.
    // This is for 3D only.
    AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM==3);

    // Forcing requires that Lx==Ly, Lz can be longer
    Real Lx = probhi[0]-problo[0];
    Real Ly = probhi[1]-problo[1];
    Real Lz = probhi[2]-problo[2];
    AMREX_ALWAYS_ASSERT(Lx==Ly);

    // Read in parameters
    ParmParse pp("turb");

    nmodes = 4;
    pp.get("nmodes", nmodes);

    div_free_force = true;
    pp.query("div_free_force", div_free_force);

    ff_factor = 4;
    pp.query("ff_factor", ff_factor);

    mode_start = 0;
    pp.query("mode_start", mode_start);

    verbose = 0;
    pp.query("verbose", verbose);

    // Inputs not yet defined. Could make runtime parameters if desired.
    int  hack_lz(0), spectrum_type(2), moderate_zero_modes(1);
    Real forcing_time_scale_min(0.5), forcing_time_scale_max(1.0), force_scale(1.0);

    // tmp CPU storage that holds everything in one flat array
    const int num_elmts=array_size*array_size*array_size;
    const int tmp_size = num_fdarray*num_elmts;
    Real tmp[tmp_size];

    // Separate out forcing data into individual Array4's
    int i_arr = 0;
    int fd_ncomp = 1;
    Dim3 fd_begin{0,0,0};
    Dim3 fd_end{array_size,array_size,array_size};

    Array4<Real> FTX(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> TAT(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPX(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPY(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPZ(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FAX(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FAY(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FAZ(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPXX(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPXY(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPXZ(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPYX(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPYY(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPYZ(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPZX(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPZY(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);
    Array4<Real> FPZZ(&tmp[(i_arr++)*num_elmts], fd_begin, fd_end, fd_ncomp);

    if (hack_lz>0) {
      if (hack_lz==1) {
        Lz = Lz/2.0;
      } else {
        Lz = Lz/hack_lz;
      }
    }

    if (verbose)
    {
        Print() << "Lx = " << Lx << std::endl;
        Print() << "Ly = " << Ly << std::endl;
        Print() << "Lz = " << Lz << std::endl;
    }

    Real Lmin = std::min(Lx,std::min(Ly,Lz));
    Real kappaMax = ((Real)nmodes)/Lmin + 1.0e-8;
    int nxmodes = nmodes*(int)(0.5+Lx/Lmin);
    int nymodes = nmodes*(int)(0.5+Ly/Lmin);
    int nzmodes = nmodes*(int)(0.5+Lz/Lmin);

    if (verbose)
    {
        Print() << "Lmin = " << Lmin << std::endl;
        Print() << "kappaMax = " << kappaMax << std::endl;
        Print() << "nxmodes = " << nxmodes << std::endl;
        Print() << "nymodes = " << nxmodes << std::endl;
        Print() << "nzmodes = " << nxmodes << std::endl;
    }

    Real freqMin = 1.0/forcing_time_scale_max;
    Real freqMax = 1.0/forcing_time_scale_min;
    Real freqDiff= freqMax-freqMin;

    if (verbose)
    {
        Print() << "forcing_time_scale_min = " << forcing_time_scale_min << std::endl;
        Print() << "forcing_time_scale_max = " << forcing_time_scale_max << std::endl;
        Print() << "freqMin = " << freqMin << std::endl;
        Print() << "freqMax = " << freqMax << std::endl;
        Print() << "freqDiff = " << freqDiff << std::endl;
    }

    // initiate the magic
    DepRand::InitRandom((unsigned long)111397);

    int mode_count = 0;

    int xstep = (int)(Lx/Lmin+0.5);
    int ystep = (int)(Ly/Lmin+0.5);
    int zstep = (int)(Lz/Lmin+0.5);

    if (verbose)
        Print() << "Mode step = " << xstep << " " << ystep << " " << zstep << std::endl;

    for (int kz = mode_start*zstep; kz <= nzmodes; kz += zstep ) {
      Real kzd = (Real)kz;
      for (int ky = mode_start*ystep; ky <= nymodes; ky += ystep ) {
        Real kyd = (Real)ky;
        for (int kx = mode_start*xstep; kx <= nxmodes; kx += xstep ) {
          Real kxd = (Real)kx;

          Real kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) );

          if (kappa<=kappaMax) {
            FTX(kx,ky,kz) = (freqMin + freqDiff*DepRand::Random() )*TwoPi;
            // Translation angles, theta=0..2Pi and phi=0..Pi
            TAT(kx,ky,kz) = DepRand::Random()*TwoPi;
            // Phases
            FPX(kx,ky,kz) = DepRand::Random()*TwoPi;
            FPY(kx,ky,kz) = DepRand::Random()*TwoPi;
            FPZ(kx,ky,kz) = DepRand::Random()*TwoPi;
            if (div_free_force==1) {
              FPXX(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPYX(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPZX(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPXY(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPYY(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPZY(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPXZ(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPYZ(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPZZ(kx,ky,kz) = DepRand::Random()*TwoPi;
            }
            // Amplitudes (alpha)
            Real thetaTmp      = DepRand::Random()*TwoPi;
            Real cosThetaTmp   = cos(thetaTmp);
            Real sinThetaTmp   = sin(thetaTmp);

            Real phiTmp        = DepRand::Random()*Pi;
            Real cosPhiTmp     = cos(phiTmp);
            Real sinPhiTmp     = sin(phiTmp);

            Real px = cosThetaTmp * sinPhiTmp;
            Real py = sinThetaTmp * sinPhiTmp;
            Real pz =               cosPhiTmp;

            Real mp2 = px*px + py*py + pz*pz;
            if (kappa < 0.000001) {
              Print() << "ZERO AMPLITUDE MODE " << kx << ky << kz << std::endl;
              FAX(kx,ky,kz) = 0.;
              FAY(kx,ky,kz) = 0.;
              FAZ(kx,ky,kz) = 0.;
            } else {
              // Count modes that contribute
              mode_count++;
              // Set amplitudes
              Real Ekh;
              if (spectrum_type==1) {
                Ekh = 1. / kappa;
              } else if (spectrum_type==2) {
                Ekh = 1. / (kappa*kappa);
              } else {
                Ekh = 1.;
              }
              if (div_free_force==1) {
                Ekh /= kappa;
              }
              if (moderate_zero_modes==1) {
                if (kx==0) Ekh /= 2.;
                if (ky==0) Ekh /= 2.;
                if (kz==0) Ekh /= 2.;
              }
              if (force_scale>0.) {
                FAX(kx,ky,kz) = force_scale * px * Ekh / mp2;
                FAY(kx,ky,kz) = force_scale * py * Ekh / mp2;
                FAZ(kx,ky,kz) = force_scale * pz * Ekh / mp2;
              } else {
                FAX(kx,ky,kz) = px * Ekh / mp2;
                FAY(kx,ky,kz) = py * Ekh / mp2;
                FAZ(kx,ky,kz) = pz * Ekh / mp2;
              }

              if (verbose)
              {
                  Print() << "Mode";
                  Print() << "kappa = " << kx << " " << ky << " " << kz << " " << kappa << " "
                          << sqrt(FAX(kx,ky,kz)*FAX(kx,ky,kz)+FAY(kx,ky,kz)*FAY(kx,ky,kz)+FAZ(kx,ky,kz)*FAZ(kx,ky,kz)) << std::endl;
                  Print() << "Amplitudes - A" << std::endl;
                  Print() << FAX(kx,ky,kz) << " " << FAY(kx,ky,kz) << " " << FAZ(kx,ky,kz) << std::endl;
                  Print() << "Frequencies" << std::endl;
                  Print() << FTX(kx,ky,kz) << std::endl;
                  Print() << "TAT" << std::endl;
                  Print() << TAT(kx,ky,kz) << std::endl;
                  Print() << "Amplitudes - AA" << std::endl;
                  Print() << FPXX(kx,ky,kz) << " " << FPYX(kx,ky,kz) << " " << FPZX(kx,ky,kz) << std::endl;
                  Print() << FPXY(kx,ky,kz) << " " << FPYY(kx,ky,kz) << " " << FPZY(kx,ky,kz) << std::endl;
                  Print() << FPXZ(kx,ky,kz) << " " << FPYZ(kx,ky,kz) << " " << FPZZ(kx,ky,kz) << std::endl;
              }
            }
          }
        }
      }
    }

    // Now let's break symmetry, have to assume high aspect ratio in z for now
    int reduced_mode_count = 0;

    for (int kz = 1; kz < zstep; kz++ ) {
      Real kzd = (Real)kz;
      for (int ky = mode_start; ky <= nymodes; ky += ystep ) {
        Real kyd = (Real)ky;
        for (int kx = mode_start; kx <= nxmodes; kx += xstep ) {
          Real kxd = (Real)kx;

          Real kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) );

          if (kappa<=kappaMax) {
            FTX(kx,ky,kz) = (freqMin + freqDiff*DepRand::Random() )*TwoPi;
            // Translation angles, theta=0..2Pi and phi=0..Pi
            TAT(kx,ky,kz) = DepRand::Random()*TwoPi;
            // Phases
            FPX(kx,ky,kz) = DepRand::Random()*TwoPi;
            FPY(kx,ky,kz) = DepRand::Random()*TwoPi;
            FPZ(kx,ky,kz) = DepRand::Random()*TwoPi;
            if (div_free_force==1) {
              FPXX(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPYX(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPZX(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPXY(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPYY(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPZY(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPXZ(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPYZ(kx,ky,kz) = DepRand::Random()*TwoPi;
              FPZZ(kx,ky,kz) = DepRand::Random()*TwoPi;
            }
            // Amplitudes (alpha)
            Real thetaTmp      = DepRand::Random()*TwoPi;
            Real cosThetaTmp   = cos(thetaTmp);
            Real sinThetaTmp   = sin(thetaTmp);

            Real phiTmp        = DepRand::Random()*Pi;
            Real cosPhiTmp     = cos(phiTmp);
            Real sinPhiTmp     = sin(phiTmp);

            Real px = cosThetaTmp * sinPhiTmp;
            Real py = sinThetaTmp * sinPhiTmp;
            Real pz =               cosPhiTmp;

            Real mp2 = px*px + py*py + pz*pz;
            if (kappa < 0.000001) {
              Print() << "ZERO AMPLITUDE MODE " << kx << ky << kz << std::endl;
              FAX(kx,ky,kz) = 0.;
              FAY(kx,ky,kz) = 0.;
              FAZ(kx,ky,kz) = 0.;
            } else {
              // Count modes that contribute
              reduced_mode_count++;
              // Set amplitudes
              Real Ekh;
              if (spectrum_type==1) {
                Ekh = 1. / kappa;
              } else if (spectrum_type==2) {
                Ekh = 1. / (kappa*kappa);
              } else {
                Ekh = 1.;
              }
              if (div_free_force==1) {
                Ekh /= kappa;
              }
              if (moderate_zero_modes==1) {
                if (kx==0) Ekh /= 2.;
                if (ky==0) Ekh /= 2.;
                if (kz==0) Ekh /= 2.;
              }
              if (force_scale>0.) {
                FAX(kx,ky,kz) = force_scale * px * Ekh / mp2;
                FAY(kx,ky,kz) = force_scale * py * Ekh / mp2;
                FAZ(kx,ky,kz) = force_scale * pz * Ekh / mp2;
              } else {
                FAX(kx,ky,kz) = px * Ekh / mp2;
                FAY(kx,ky,kz) = py * Ekh / mp2;
                FAZ(kx,ky,kz) = pz * Ekh / mp2;
              }

              if (verbose)
              {
                  Print() << "Mode";
                  Print() << "kappa = " << kx << " " << ky << " " << kz << " " << kappa << " "
                          << sqrt(FAX(kx,ky,kz)*FAX(kx,ky,kz)+FAY(kx,ky,kz)*FAY(kx,ky,kz)+FAZ(kx,ky,kz)*FAZ(kx,ky,kz)) << std::endl;
                  Print() << "Amplitudes - A" << std::endl;
                  Print() << FAX(kx,ky,kz) << " " << FAY(kx,ky,kz) << " " << FAZ(kx,ky,kz) << std::endl;
                  Print() << "Frequencies" << std::endl;
                  Print() << FTX(kx,ky,kz) << std::endl;
                  Print() << "TAT" << std::endl;
                  Print() << TAT(kx,ky,kz) << std::endl;
                  Print() << "Amplitudes - AA" << std::endl;
                  Print() << FPXX(kx,ky,kz) << " " << FPYX(kx,ky,kz) << " " << FPZX(kx,ky,kz) << std::endl;
                  Print() << FPXY(kx,ky,kz) << " " << FPYY(kx,ky,kz) << " " << FPZY(kx,ky,kz) << std::endl;
                  Print() << FPXZ(kx,ky,kz) << " " << FPYZ(kx,ky,kz) << " " << FPZZ(kx,ky,kz) << std::endl;
              }
            }
          }
        }
      }
    }

    Print() << "mode_count = " << mode_count << std::endl;
    Print() << "reduced_mode_count = " << reduced_mode_count << std::endl;
    if (spectrum_type==1) {
      Print() << "Spectrum type 1" << std::endl;
    } else if (spectrum_type==2) {
      Print() << "Spectrum type 2" << std::endl;
    } else {
      Print() << "Spectrum type OTHER" << std::endl;
    }

// Now allocate forcedata and copy in tmp array.
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        forcedata = static_cast<Real*>(The_Arena()->alloc(tmp_size*sizeof(Real)));
        Gpu::htod_memcpy_async(forcedata, tmp, tmp_size*sizeof(Real));
    }
    else
#endif
    {
        forcedata = static_cast<Real*>(The_Pinned_Arena()->alloc(tmp_size*sizeof(Real)));
        std::memcpy(forcedata, tmp, tmp_size*sizeof(Real));
    }
}
