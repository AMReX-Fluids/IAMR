..
  This needs to be reorganized, and probably some parts should be moved to different Chapters


Code structure
==============

The code structure in the IAMR/ directory is as follows:

-  Source/: source code

-  Exec/: various problem run directories, including:

   -  run2d/

   -  run3d/

   -  eb_run2d/

   -  eb_run3d/

   -  run_2d_particles/

-  Docs/: you’re reading this now!

An Overview of IAMR
===================

IAMR is built upon the AMReX C++ framework. This provides
high-level classes for managing an adaptive mesh refinement simulation,
including the core data structures required in AMR calculations.
Since IAMR leaverages heavily from the AMReX library,
it’s documentation :ref:`amrex:amrex_doc_indx` 
is a useful resource in addition to this User’s Guide.

..
    (at https://amrex-codes.github.io/amrex/docs_html/index.html)


The IAMR simulation begins in IAMR/Source/main.cpp where an instance
of the AMReX Amr class is created:

::

      Amr* amrptr = new Amr;

The initialization, including calling a problem’s initdata()
routine and refining the base grid occurs next through

::

      amrptr->init(strt_time,stop_time);

And then comes the main loop over coarse timesteps until the
desired simulation time is reached:

::

      while ( amrptr->okToContinue()                            &&
             (amrptr->levelSteps(0) < max_step || max_step < 0) &&
             (amrptr->cumTime() < stop_time || stop_time < 0.0) )

      {
         //
         // Do a timestep.
         //
         amrptr->coarseTimeStep(stop_time);
      }

This uses the AMReX machinery to do the necessary subcycling in time,
including synchronization between levels, to advance the level hierarchy
forward in time.

State Data
~~~~~~~~~~

The StateData class structure defined by AMReX is the data container
used to store the field data associated with the state on a single AMR level
during an IAMR run. The entire state consists of a dynamic union, or hierarchy, of
nested StateData objects. Periodic regrid operations modify the hierarchy,
changing the shape of the data containers at the various levels according to
user-specified criteria; new StateData objects are created
for the affected levels, and are filled with the “best” (finest) available
data at each location. Instructions for building and managing StateData are
encapsulated in the AMReX class, StateDescriptor; as discussed later,
a StateDescriptor will be created for each type of state field, and
will include information about data centering, required grow cells, and
instructions for transferring data between AMR levels during various synchronization
operations.

In IAMR/Source/NavieStokesBase.H, the enum StateType defines the
different state descriptors for IAMR. These are setup during the
run by code in NS_setup.cpp, and include (but are not limited to):

-  State_Type: the cell-centered density, velocity, and other scalars (tracers)

-  Press_Type: the node-centered dynamic pressure field.

-  Divu_Type: Stores the right-hand-side of the constraint
   (only matters for low Mach flows when this is nonzero).

-  Dsdt_Type: Stores the time-derivative of the right-hand-side of the constraint
   (only matters for low Mach flows when this is nonzero).

Each StateData object has two MultiFabs, one each for
old and new times, and can provide an interpolated copy of the state at any time between the two.
Alternatively, can also access the data containers directly, for instance:

::

    MultiFab& S_new = get_new_data(State_Type);

gets a pointer to the multifab containing the hydrodynamics state data
at the new time (here State_Type is the enum defined in
NavierStokesBase.H) (note that the class NavierStokes
is a derived classes of NavierStokesBase).

MultiFab data is distributed in space at the granularity of
each Box in its BoxArray. We iterate over MultiFabs using a special
iterator, MFIter, which knows about the locality of the data—only the boxes owned by the
processor will be included in the loop on each processor. An example loop
(taken from code in NavierStokesBase.cpp):

::

    //
    // Fill rho at half time
    //
    #ifdef _OPENMP
    #pragma omp parallel if (Gpu::notInLaunchRegion())
    #endif
        for (MFIter mfi(rho_half,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox();

            // half time
            auto const& rho_h = rho_half.array(mfi);
            // previous time
            auto const& rho_p = rho_ptime.array(mfi);
            // current time
            auto const& rho_c = rho_ctime.array(mfi);

            amrex::ParallelFor(bx, [rho_h, rho_p, rho_c] 
            AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
               rho_h(i,j,k) = 0.5 * (rho_p(i,j,k) + rho_c(i,j,k));
            });
        }
    }

Here, ++mfi iterates to the next FArrayBox owned by the MultiFab,
and mfi.isValid() returns false after we’ve reached
the last box contained in the MultiFab, terminating the loop.
rho_half.array(mfi) creates an object for accessing FArrayBox data in
a more array like manner using operator()
(more details are in the AMReX documentation at
https://amrex-codes.github.io/amrex/docs_html/Basics.html#sec-basics-array4 ).

Here ParallelFor takes two arguments. The first argument is a Box specifying the iteration index space, and the second argument is a C++ lambda function that works on cell (i,j,k). Variables rho_half, rho_ptime and rho_ctime in the lambda function are captured by value from the enclosing scope. The code above is performance portable. It works with and without GPU support. When IAMR is built with GPU support (USE_CUDA=TRUE), AMREX_GPU_DEVICE indicates that the lambda function is a device function and ParallelFor launches a GPU kernel to do the work. When it is built without GPU support, AMREX_GPU_DEVICE has no effects whatsoever. It should be emphasized that ParallelFor does not start an OpenMP parallel region. The OpenMP parallel region will be started by the pragma above the MFIter loop if it is built with OpenMP and without enabling GPU (USE_OMP=TRUE and USE_CUDA=TRUE are not compatible). Tiling is turned off if GPU is enabled so that more parallelism is exposed to GPU kernels. Also note that when tiling is off, tilebox returns validbox.
(more details are in the AMReX documentation at
https://amrex-codes.github.io/amrex/docs_html/GPU.html#sec-gpu-for ).

Setting Up Your Own Problem
===========================

To define a new problem, we create a new inputs file in
a run directory and modify
IAMR/Source/prob/prob_init.cpp accordingly.
The simplest way to get started is to copy the inputs files from an existing
problem. Here we describe how to customize your problem.

There are several files involved in setting up an IAMR problem. It’s possible to
create your own new setup by only touching the first of these (prob_initData())
and changing parameters through the inputs file (see section [sec:inputs]).
Here we list the most relevant problem
setup files and thier purpose. If you need further help setting up your problem, please
contact us.

-  prob_initData():
   Read in initial conditions and problem parameters from the inputs file,
   and initialize the state data (velocity, density, etc.).

-  NS_error.cpp: Define the error estimation criteria used for tagging cells for
   refinement.
   More details in section [sec:tagging]

-  NS_setup.cpp: Declare state and derived variables.
   Specify how to fill grow cells for each state or derived variable.
   More details in sections [sec:boundaries]

-  NS_derive.cpp: Define derived variables.
   More details in sections [sec:derivedVariables]

-  NS_BC.H: Define the mapping from physical boundary conditions (e.g. outflow)
   to mathematical (e.g. first order extrapolation from last interior cell).
   More details in section [sec:physicalBCs]

-  NS_bcfill.H:
   Define the boundary filling functions for external Dirichlet (i.e. user supplied)
   boundary conditions. Constant Dirichlet conditions can be specified in the
   inputs file without needing to alter NS_bcfill.H.
   More details in section [sec:physicalBCs]

Boundaries
==========

In AMReX, we are primarily concerned with enabling structured-grid
computations. A key aspect of this is the use of “grow” cells
around the “valid box” of cells over which we wish to apply stencil operations.
Grow cells, filled properly, are conveniently located temporary
data containers that allow us to separate the steps of data preparation
(including communication, interpolation, or other complex manipulation)
from stencil application. The steps that are required to fill grow cells
depends on where the cells “live” in the computational domain.

Boundaries Between Grids and Levels
-----------------------------------

Most of our state data is cell-centered, and often the grow cells are
as well. When the cells lie directly over cells of a neighboring box
at the same AMR refinement level, these are “fine-fine” cells, and are
filled by direct copy (including any MPI communication necessary to enable
that copy). Note that fine-fine boundary also include grow cells that
cover valid fine cells through a periodic boundary.

When the boundary between valid and grow cells is coincident
with a coarse-fine boundary, these coarse-fine grow cells will hold cell-centered
temporary data that generated by interpolation (in space and time) of the
underlying coarse data. This operation requires auxiliary metadata to define
how the interpolation is to be done, in both space and time. Importantly,
the interpolation also requires that coarse data be well-defined over
a time interval that brackets the time instant for which we are evaluating
the grow cell value – this places requirements on how the time-integration
of the various AMR levels are sequenced relative to eachother.
In AMReX, the field data associated with the system state, as well as the metadata
associated with inter-level transfers, is bundled (encapsulated) in
a class called “StateData”. The metadata
is defined in NS_setup.cpp – search for
cell_cons_interp, for example – which is “cell conservative
interpolation”, i.e., the data is cell-based (as opposed to node-based
or edge-based) and the interpolation is such that the average of the
fine values created is equal to the coarse value from which they came.
(This wouldn’t be the case with straight linear interpolation, for
example.) A number of interpolators are provided with AMReX and
user-customizable ones can be added on the fly.


.. _sec:physicalBCs:

Physical Boundaries
-------------------

The last type of grow cell exists at physical boundaries. These are special for
a couple of reasons. First, the user must explicitly specify how they are to be
filled, consistent with the problem being run. AMReX provides a number of
standard condition types typical of PDE problems (reflecting, extrapolated, etc),
and a special one that indicates external Dirichlet. In the case of Dirichlet,
the user supplies data to fill grow cells.

IAMR provides the ability to specify constant Dirichlet BCs
in the inputs file (see section [sec:dirichlet]).
Users can create more complex Dirichlet boundary condtions by writing
their own fill function in NS_bcfill.H, then using that function to create
an amrex::StateDescriptor::BndryFunc object and specifying which variables
will use it in NS_setup.cpp.

It is important to note that external Dirichlet boundary data is to be specified as
if applied on the face of the cell bounding the domain, even for cell-centered
state data. For cell-centered data, the array passed into the
boundary condition code is filled with cell-centered values in the valid
region and in fine-fine, and coarse-fine grow cells. Additionally, grow cells
for standard extrapolation and reflecting boundaries are pre-filled. The
differential operators throughout IAMR are aware of the special boundaries
that are Dirichlet and wall-centered, and the stencils are adjusted accordingly.

For convenience, IAMR provides a limited set of mappings from a physics-based boundary condition
specification to a mathematical one that the code can apply. This set can be extended
by adjusting the corresponding translations in NS_BC.H, but, by default, includes
(See AMReX/Src/Base/AMReX_BC_TYPES.H for more detail):

-  *Outflow*:

   -  velocity: FOEXTRAP

   -  temperature: FOEXTRAP

   -  scalars: FOEXTRAP

-  *No Slip Wall with Adiabatic Temp*:

   -  velocity: EXT_DIR, :math:`u=v=0`

   -  temperature: REFLECT_EVEN, :math:`dT/dt=0`

   -  scalars: HOEXTRAP

-  *Slip Wall with Adiabatic Temp*:

   -  velocity: EXT_DIR, :math:`u_n=0`; HOEXTRAP, :math:`u_t`

   -  temperature: REFLECT_EVEN, :math:`dT/dn=0`

   -  scalars: HOEXTRAP

The keywords used above are defined:

-  INT_DIR: data taken from other grids or interpolated

-  EXT_DIR: data specified on EDGE (FACE) of bndry

-  HOEXTRAP: higher order extrapolation to EDGE of bndry

-  FOEXTRAP: first order extrapolation from last cell in interior

-  REFLECT_EVEN: :math:`F(-n) = F(n)` true reflection from interior cells

-  REFLECT_ODD: :math:`F(-n) = -F(n)` true reflection from interior cells

Derived Variables
=================

IAMR has the ability to created new variables derived from the state variables.
A few derived variables are provided with IAMR, which can be used as examples for
creating user defined derived variables.
Users create derived variables by adding a function to create them in
NS_derive.H and NS_derive.cpp, and then adding the variable to the
derive_lst in NS_setup.cpp.

Access to the derived variable is through one of two amrex:AmrLevel functions
(which are inherited by NavierStokesBase and NavierStokes):

::

        /**
        * \brief Returns a MultiFab containing the derived data for this level.
        * The user is responsible for deleting this pointer when done
        * with it.  If ngrow>0 the MultiFab is built on the appropriately
        * grown BoxArray.
        */
        virtual std::unique_ptr<MultiFab> derive (const std::string& name,
                              Real               time,
                              int                ngrow);
        /**
        * \brief This version of derive() fills the dcomp'th component of mf
        * with the derived quantity.
        */
        virtual void derive (const std::string& name,
                             Real               time,
                             MultiFab&          mf,
                             int                dcomp);

As an example, mag\_vort is a derived variable provided with IAMR, which
returns the magnitude of the vorticity of the flow.
A multifab filled with the magnitude of the vorticity can be obtained via

::

      std::unique_ptr<MultiFab> vort;
      vort = derive(mag_vort, time, ngrow);
      //
      // do something with vorticity...
      //
      vort.reset();  

The FillPatchIterator
---------------------

A FillPatchIterator is a AMReX object tasked with the job of
filling rectangular patches of state data, possibly including grow cells,
and, if so, utilizing all the metadata discussed above that is provided by
the user. Thus, a FillPatchIterator can only be constructed on
a fully registered StateData object, and is the preferred
process for filling grown platters of data prior to most stencil
operations (e.g., explicit advection operators, which may require
several grow cells). It should be mentioned that a FillPatchIterator
fills temporary data via copy operations, and therefore does not
directly modify the underlying state data. In the code, if the state
is modified (e.g., via an advective “time advance”, the new data
must be copied explicitly back into the StateData containers.

Use of FillPatchIterator as an iterator has been depreciated in favor
of MFIter, which supports tiling (see section :ref:`Chap:Parallel`).
However, IAMR continues to use
FillPatchIterator for creating temporaries with filled grow cells.

For example, the following code demonstrates the calling sequence to
create and use a FillPatchIterator for preparing a rectangular patch of
data that includes the “valid region” plus NUM_GROW grow cells. Here,
the valid region is specified as a union of rectangular boxes making up the
box array underlying the MultiFab S_new, and NUM_GROW cells are
added to each box in all directions to create the temporary patches to
be filled.

::

      FillPatchIterator fpi(*this, S_new, NUM_GROW,
                            time, State_Type, strtComp, NUM_STATE);
      // Get a reference to the temporary platter of grown data
      MultiFab& S = fpi.get_mf();

Here the FillPatchIterator fills the patch
with data of type “State_Type” at time “time”,
starting with component strtComp and including a total of
NUM_STATE components. When the FillPatchIterator goes out of scope, it
and the temporary data platters are destroyed (though much of the
metadata generated during the operation is cached internally
for performance). Notice that since NUM_GROW can be any
positive integer (i.e., that the grow region can extend over an arbitrary
number of successively coarser AMR levels), this key operation can hide an
enormous amount of code and algorithm complexity.

Parallel I/O
============

Both checkpoint files and plotfiles are actually folders containing
subfolders: one subfolder for each level of the AMR hierarchy.
The fundamental data structure we read/write to disk is a MultiFab,
which is made up of multiple FAB’s, one FAB per grid. Multiple
MultiFabs may be written to each folder in a checkpoint file.
MultiFabs of course are shared across CPUs; a single MultiFab may be
shared across thousands of CPUs. Each CPU writes the part of the
MultiFab that it owns to disk, but they don’t each write to their own
distinct file. Instead each MultiFab is written to a runtime
configurable number of files N (N can be set in the inputs file as the
parameter amr.checkpoint_nfiles and amr.plot_nfiles; the
default is 64). That is to say, each MultiFab is written to disk
across at most N files, plus a small amount of data that gets written
to a header file describing how the file is laid out in those N files.

What happens is :math:`N` CPUs each opens a unique one of the :math:`N` files into
which the MultiFab is being written, seeks to the end, and writes
their data. The other CPUs are waiting at a barrier for those :math:`N`
writing CPUs to finish. This repeats for another :math:`N` CPUs until all the
data in the MultiFab is written to disk. All CPUs then pass some data
to CPU 0 which writes a header file describing how the MultiFab is
laid out on disk.

We also read MultiFabs from disk in a “chunky” manner, opening only :math:`N`
files for reading at a time. The number :math:`N`, when the MultiFabs were
written, does not have to match the number :math:`N` when the MultiFabs are
being read from disk. Nor does the number of CPUs running while
reading in the MultiFab need to match the number of CPUs running when
the MultiFab was written to disk.

Think of the number :math:`N` as the number of independent I/O pathways in
your underlying parallel filesystem. Of course a “real” parallel
filesytem should be able to handle any reasonable value of :math:`N`. The
value -1 forces :math:`N` to the number of CPUs on which you’re
running, which means that each CPU writes to a unique file, which can
create a very large number of files, which can lead to inode issues.
