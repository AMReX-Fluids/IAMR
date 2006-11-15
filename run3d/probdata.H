      REAL_T  denerr, vorterr, adverr, temperr, randfact, bubgrad
      REAL_T  denfact, xblob, yblob, zblob, radblob, velfact 
      REAL_T  rhozero, tempzero, r_d, c_d, adv_vel, radvort
      REAL_T  den1,den2,vel1,vel2,delta0,xlev1,zlev1,amag,freq(10),mag(10)
      REAL_T  f_problo(SDIM), f_probhi(SDIM)
      REAL_T  turb_scale,domnlo(SDIM),domnhi(SDIM)
      REAL_T  grav_angle, omega
      REAL_T  AXY(10), BXY(10), CXY(10), PXY(10), QXY(10), RXY(10)
      REAL_T  AZX(10), BZX(10), CZX(10), PZX(10), QZX(10), RZX(10)
      REAL_T  AYZ(10), BYZ(10), CYZ(10), PYZ(10), QYZ(10), RYZ(10)
      integer probtype, adv_dir, axis_dir, vb_unifdir, nmodes
      common /probdt1/ probtype, adv_dir, axis_dir, nmodes
      common /probdt2/  den1,den2,vel1,vel2,delta0,
     &                  xlev1,zlev1,freq,mag, amag
      common /probdt3/ vorterr, adverr, denerr, temperr
      common /probdt4/ denfact, xblob, yblob, zblob, radvort,
     &                 radblob, velfact, randfact, bubgrad
      common /probdt5/ rhozero, tempzero, r_d, c_d, adv_vel
      common /probdt6/ f_problo, f_probhi,domnlo,domnhi,turb_scale
      common /probdt7/ vb_unifdir
      common /probdt9/ grav_angle, omega,
     &                 AXY, BXY, CXY, PXY, QXY, RXY,
     &                 AZX, BZX, CZX, PZX, QZX, RZX,
     &                 AYZ, BYZ, CYZ, PYZ, QYZ, RYZ

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
