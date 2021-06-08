.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _sec:EB-basics:

Constructing Embedded Boundaries in incflo
============================================

MincfloFiX uses AMReX's constructive solid geometry framework defined in the namespace
:cpp:`amrex::EB2`. See the `AMReX EB documentation`_ for more details. These are
defined in ``incflo/src/embedded_boundaries ``. 

How incflo Constructs the EB Geometry
---------------------------------------

Once a geometry is selected by :cpp:`incflo::make_eb_geometry`, the procedure is
the same for (almost) all geometries. Also see the `AMReX geometry
documentation`_ for information on how to construct new geometries:

1. Construct an implicit function representing the geometry (using the language
   of constructive solid geometry). For example

.. highlight:: c++

::

   EB2::CylinderIF my_cyl(radius, height, direction, center, inside);
   auto gshop_cyl = EB2::makeShop(my_cyl);

2. Call :cpp:`incflo::build_eb_levels(gshop)` this function builds the EB levels
   and fills the implicit function :cpp:`MultiFab` (the later being used to
   construct the level-set function). 

incflo's EB Data Structures
-------------------------

The :cpp:`incflo` class stores the following EB data:

.. highlight:: c++

::

   //! EB levels representing fluid boundary conditions
   Vector<const EB2::Level *> eb_levels;

   //! EB factory that lives on the fluid grids
   Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > ebfactory;

A note about constructing EB Levels
-----------------------------------

incflo builds EB levels in :cpp:`incflo::build_eb_levels` (via
:cpp:`LSCore<F>::BuildEBLevel`)

.. highlight:: c++

::

   EB2::Build(gshop, geom[lev], required_crse_lev, max_crse_level);
   const EB2::IndexSpace & ebis = EB2::IndexSpace::top();


When building an EB level, the maximum coarsening level (:cpp:`int
max_crse_level`) and the required coarsening level (:cpp:`int
required_crse_lev`) need to be specified. The reason for this is that we need to
specify to which level of coarseness the EB is still defined. It might not be
immediately obvious, but the Poisson solver (used in the fluid solve) also
depends indirectly on these parameters. Thus changing these during EB level
creation might restrict how many levels the MLMG solver can use, and therefore
give slightly different answers in the fluid solve.

incflo Initialization Process
-------------------------------

Since incflo requires the volume fraction when building grids (because this is
needed by :cpp:`incflo::ErrorEst`), the EB geometries need to be built before
calling :cpp:`incflo::Init`. The recommended procedure therefore is

.. highlight:: c++

::

   // Default constructor (geom[lev] is defined here)
   incflo my_incflo;

   // Initialize internals from ParamParse database
   my_incflo.InitParams(solve_fluid, solve_dem, call_udf);

   // Initialize memory for data-array internals
   my_incflo.ResizeArrays();

   // Construct EB (must be done _before_ incflo::Init)
   my_incflo.make_eb_geometry();

   // Initialize derived internals. Grids are create here.
   my_incflo.Init(dt, time);

   // Create EB factories on new grids
   my_incflo.make_eb_factories();

   // Finish constructing levels
   my_incflo.InitLevelData(dt,time);

   // Regrid (ensure all MultiFabs are on their correct grids)
   my_incflo.Regrid();

The grids for each level are build in the :cpp:`incflo::Init` by invoking the
initialization functions inherited from :cpp:`amrex::AmrCore`.

.. highlight:: c++

::

   // This tells the AmrMesh class not to iterate when creating the initial
   // grid hierarchy
   SetIterateToFalse();

   // This tells the Cluster routine to use the new chopping routine which
   // rejects cuts if they don't improve the efficiency
   SetUseNewChop();

   // This Builds the new Grids
   InitFromScratch(0.);
