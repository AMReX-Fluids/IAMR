The inputs.3d.multi_chamber case sets up a problem with a verical pipe from the YHI boundary that dumps to a periodic chamber that has a drain.  The drain dumps to a lower chamber in contact with an outflow boundary on YLO.  Currently, this case fails at step 857, where the level project diverges (after running along just fine up to this point).  The inputs can be extended to make the lower chamber very long in order to avoid backflow on the outflow face and the code still crashes in the same way.

This setup comes from Sigfried Haering at SNL.
-MD, 10/30/2020

