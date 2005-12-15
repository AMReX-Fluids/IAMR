      REAL_T  denerr, vorterr, adverr, temperr, randfact, bubgrad
      REAL_T  denfact, xblob, yblob, zblob, radblob, velfact 
      REAL_T  rhozero, tempzero, r_d, c_d, adv_vel, radvort
      REAL_T  den1,den2,vel1,vel2,delta0,xlev1,zlev1,amag,freq(10),mag(10)
      REAL_T f_problo(SDIM), f_probhi(SDIM), turb_scale
      integer probtype, adv_dir, axis_dir
      common /probdt1/ probtype, adv_dir, axis_dir
      common /probdt2/  den1,den2,vel1,vel2,delta0,
     &                  xlev1,zlev1,freq,mag
      common /probdt3/ vorterr, adverr, denerr, temperr
      common /probdt4/ denfact, xblob, yblob, zblob, radvort,
     &                 radblob, velfact, randfact, bubgrad
      common /probdt5/ rhozero, tempzero, r_d, c_d, adv_vel
      common /probdt6/f_problo, f_probhi, turb_scale

      logical forceInflow
      REAL_T probSizeFile(3)
      common /inflow1/forceInflow
      common /inflow2/probSizeFile
      save /inflow1/, /inflow2/

