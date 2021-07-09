.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _sec:EB-basics:

Constructing Embedded Boundaries in IAMR
============================================

IAMR uses AMReX's constructive solid geometry framework defined in the namespace
:cpp:`amrex::EB2`. See the AMReX EB documentation (:ref:`amrex:sec:EB:ebinit`) for more details.

How IAMR Constructs the EB Geometry
---------------------------------------

Once a geometry is selected by :cpp:`IAMR::make_eb_geometry`, the procedure is
the same for (almost) all geometries. Also see the AMReX geometry
documentation :ref:`amrex:sec:EB:ebinit:IF` for information on how to construct new geometries:

1. Construct an implicit function representing the geometry (using the language
   of constructive solid geometry). For example

.. highlight:: c++

::

   EB2::CylinderIF my_cyl(radius, height, direction, center, inside);
   auto gshop_cyl = EB2::makeShop(my_cyl);

2. Call :cpp:`IAMR::build_eb_levels(gshop)` this function builds the EB levels
   and fills the implicit function :cpp:`MultiFab` (the later being used to
   construct the level-set function). 

