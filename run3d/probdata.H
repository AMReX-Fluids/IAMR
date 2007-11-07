      REAL_T  denerr, vorterr, adverr, temperr, randfact, bubgrad
      REAL_T  denfact, xblob, yblob, zblob, radblob, velfact 
      REAL_T  rhozero, tempzero, r_d, c_d, adv_vel, radvort
      REAL_T  den1,den2,vel1,vel2,delta0,xlev1,zlev1,amag,freq(10),mag(10)
      REAL_T  f_problo(SDIM), f_probhi(SDIM)
      REAL_T  turb_scale,domnlo(SDIM),domnhi(SDIM)
      REAL_T  grav_angle, omega, infl_time_offset, injection_time
      REAL_T  jet_x, jet_y, jet_width, jet_rho, jet_vel, coflow_rho, coflow_vel, jet_temp
      REAL_T  ref_centre, ref_radius, ref_height
      REAL_T  thermal_expansion, heating_coeff
      integer probtype, adv_dir, axis_dir, vb_unifdir, blrandseed
      integer ref_type
      common /probdt1/ probtype, adv_dir, axis_dir, ref_type
      common /probdt2/ den1,den2,vel1,vel2,delta0,
     &                 xlev1,zlev1,freq,mag, amag, injection_time
      common /probdt3/ vorterr, adverr, denerr, temperr, ref_centre, ref_radius
      common /probdt4/ denfact, xblob, yblob, zblob, radvort,
     &                 radblob, velfact, randfact, bubgrad
      common /probdt5/ rhozero, tempzero, r_d, c_d, adv_vel
      common /probdt6/ f_problo, f_probhi,domnlo,domnhi,turb_scale
      common /probdt7/ vb_unifdir, blrandseed
      common /probdt9/ grav_angle, omega, infl_time_offset
      common /probdt10/ jet_x, jet_y, jet_width, jet_rho, jet_vel, coflow_rho, coflow_vel, jet_temp, ref_height
      common /probdt11/ thermal_expansion, heating_coeff

      integer rt_max_freq
      parameter (RT_MAX_FREQ = 10)	
      REAL_T rt_splitx, rt_xfrontw, rt_den_1, rt_den_2
      REAL_T rt_ranampl(RT_MAX_FREQ,RT_MAX_FREQ)
      REAL_T rt_ranphse(RT_MAX_FREQ,RT_MAX_FREQ,2)
      REAL_T rt_pertamp
      REAL_T rt_graddenerr
      integer rt_nfreq
      common /probdt8/ rt_splitx, rt_xfrontw, rt_den_1, rt_den_2
      common /probdt8/ rt_ranampl, rt_ranphse, rt_pertamp
      common /probdt8_i/ rt_nfreq
      common /probdt8/ rt_graddenerr

      REAL_T rt_PI
      parameter (rt_PI = 3.1415926535897932385)
      logical forceInflow
      REAL_T probSizeFile(3)
      common /inflow1/forceInflow
      common /inflow2/probSizeFile
      save /inflow1/, /inflow2/

      REAL_T Vco, Rfu, Rtran, tVco_l, tVco_r, Vco_l, Vco_r
      common /inflow3/Vco, Rfu, Rtran, tVco_l, tVco_r, Vco_l, Vco_r
      save /inflow3/
