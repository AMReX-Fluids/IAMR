c
c ::: Inflow Forcing Common Blocks :::
c
      integer strmwse_dir, numInflPlanesStore
      DOUBLE PRECISION  convVel
      logical forceLo, forceHi
      character flct_file*(80)

      common /InflFrc_Comms/convVel, strmwse_dir, numInflPlanesStore,
     $                      forceLo, forceHi
      common /InflFrc_File/ flct_file
      save /InflFrc_Comms/, /InflFrc_File/

