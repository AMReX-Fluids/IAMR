      REAL_T  forcing_xlength, forcing_ylength, forcing_zlength, forcing_epsilon
      REAL_T  forcing_time_scale_min, forcing_time_scale_max, force_scale
      REAL_T  AXY(32), BXY(32), CXY(32), DXY(32), PXY(32), QXY(32), RXY(32)
      REAL_T  AZX(32), BZX(32), CZX(32), DZX(32), PZX(32), QZX(32), RZX(32)
      REAL_T  AYZ(32), BYZ(32), CYZ(32), DYZ(32), PYZ(32), QYZ(32), RYZ(32)
      REAL_T  FTX(0:32,0:32,0:32), FTY(0:32,0:32,0:32), FTZ(0:32,0:32,0:32), TAT(0:32,0:32,0:32), TAP(0:32,0:32,0:32)
      REAL_T  FPX(0:32,0:32,0:32), FPY(0:32,0:32,0:32), FPZ(0:32,0:32,0:32), FAX(0:32,0:32,0:32), FAY(0:32,0:32,0:32), FAZ(0:32,0:32,0:32)
      REAL_T  FPXX(0:32,0:32,0:32), FPYX(0:32,0:32,0:32), FPZX(0:32,0:32,0:32)
      REAL_T  FPXY(0:32,0:32,0:32), FPYY(0:32,0:32,0:32), FPZY(0:32,0:32,0:32)
      REAL_T  FPXZ(0:32,0:32,0:32), FPYZ(0:32,0:32,0:32), FPZZ(0:32,0:32,0:32)
      integer nmodes, nxmodes, nymodes, nzmodes, mode_start, hack_lz
      integer spectrum_type, forcing_type, use_rho_in_forcing, do_mode_division, div_free_force, forcing_twice_wavelength, moderate_zero_modes
      common /forcedti/nmodes, nxmodes, nymodes, nzmodes, mode_start, hack_lz, &
                      spectrum_type, forcing_type,use_rho_in_forcing, do_mode_division, div_free_force, forcing_twice_wavelength, moderate_zero_modes
      common /forcedt/ force_scale, forcing_epsilon,  &
                      forcing_xlength, forcing_ylength, forcing_zlength,&
                      forcing_time_scale_min, forcing_time_scale_max,&
                      AXY, BXY, CXY, DXY, PXY, QXY, RXY,&
                      AZX, BZX, CZX, DZX, PZX, QZX, RZX,&
                      AYZ, BYZ, CYZ, DYZ, PYZ, QYZ, RYZ,&
                      FTX, FTY, FTZ, TAT, TAP,&
                      FPX, FPY, FPZ, FAX, FAY, FAZ,&
                      FPXX, FPYX, FPZX, FPXY, FPYY, FPZY, FPXZ, FPYZ, FPZZ
