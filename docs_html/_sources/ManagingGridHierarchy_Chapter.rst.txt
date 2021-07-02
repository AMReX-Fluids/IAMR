.. role:: cpp(code)
   :language: c++

Gridding and Load Balancing
===========================

IAMR has a great deal of flexibility when it comes to how to decompose the
computational domain into individual rectangular grids, and how to distribute
those grids to MPI ranks.  There can be grids of different sizes, 
more than one grid per MPI rank, and different strategies for distributing the grids to MPI ranks.
IAMR relies on AMReX for the implementation. For more information, please see AMReX's documentation,
found here: :ref:`amrex:Chap:ManagingGridHierarchy`.

See :ref:`sec:gridCreation` and :ref:`Chap:InputsLoadBalancing` for how grids are created, i.e. how the :cpp:`BoxArray` on which 
:cpp:`MultiFabs` will be built is defined at each level.

See :ref:`amrex:sec:load_balancing` for the strategies AMReX supports for distributing
grids to MPI ranks, i.e. defining the :cpp:`DistributionMapping` with which 
:cpp:`MultiFabs` at that level will be built.  

When running on multicore machines with OpenMP, we can also control the distribution of 
work by setting the size of grid tiles (by defining :cpp:`fabarray.mfiter_tile_size`). 
We can also specify the strategy for assigning tiles to OpenMP threads.  
See :ref:`Chap:Parallel` and AMReX's docmumentation (:ref:`amrex:sec:basics:mfiter:tiling`) for more about tiling.

