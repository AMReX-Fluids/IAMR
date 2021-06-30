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

