.. role:: cpp(code)
   :language: c++

.. _Chap:parallel:
	      
Parallelization
===============

AMReX uses a hybrid MPI + X approach to parallelism,
where X = OpenMP for multicore machines, and CUDA/HIP/DCP++ for CPU/GPU systems.
The basic idea is that MPI is used to distribute individual boxes across
nodes while X is used to distribute the work in local boxes
within a node. The OpenMP approach in AMReX is optionally
based on *tiling* the box-based data structures. Both the tiling and
non-tiling approaches to work distribution are discussed below. Also see 
the discussion of tiling in AMReX's documentation, :ref:`amrex:sec:basics:mfiter`.


AMReX’s Non-Tiling Approach
---------------------------

At the highest abstraction level, we have MultiFab (mulitple
FArrayBoxes). A MultiFab contains an array of Boxes (a Box contains integers specifying the index space it
covers), including Boxes owned by other processors for the
purpose of communication, an array of MPI ranks specifying which MPI
processor owns each Box, and an array of pointers to FArrayBoxes owned by this MPI processor.
A typical usage of MultiFab is as follows,

::

      for (MFIter mfi(mf); mfi.isValid(); ++mfi) // Loop over boxes
      {
        // Get the index space of this iteration
        const Box& box = mfi.validbox(); 

        // Get a reference to mf's data as a multidimensional array 
        auto& data = mf.array(mfi);

        // Loop over "box" to update data.
        // On CPU/GPU systems, this loop executes on the GPU
      }

A few comments about this code

-  Here the iterator, :cpp:`mfi`, will perform the loop only over the
   boxes that are local to the MPI task. If there are 3 boxes on the
   processor, then this loop has 3 iterations.

-  box as returned from :cpp:`mfi.validbox()` does not include
   ghost cells. We can get the indices of the valid zones as box.loVect and box.hiVect.

AMReX’s Tiling Approach
-----------------------

There are two types of tiling that people discuss. In *logical
tiling*, the data storage in memory is unchanged from how we do things
now in pure MPI. In a given box, the data region is stored
contiguously). But when we loop in OpenMP over a box, the tiling
changes how we loop over the data. The alternative is called *separate tiling*—here the data storage in memory itself is changed
to reflect how the tiling will be performed. This is not considered
in AMReX.

In our logical tiling approach, a box is logically split into tiles,
and a MFIter loops over each tile in each box. Note that the
non-tiling iteration approach can be considered as a special case of
tiling with the tile size equal to the box size.

::

      bool tiling = true;
      for (MFIter mfi(mf,tiling); mfi.isValid(); ++mfi) // Loop over tiles
      {
        // Get the index space of this iteration
        const Box& box = mfi.tilebox(); 

        // Get a reference to mf's data as a multidimensional array 
        auto& data = mf.array(mfi);

        // Loop over "box" to update data.
      }

Note that the code is almost identical to the one with the non-tiling approach.
Some comments:

-  The iterator now takes an extra argument to turn on tiling
   (set to true). There is another interface fo MFIter
   that can take an IntVect that explicitly gives the tile size
   in each coordinate direction.

   If we don’t explictly specify the tile size at the loop, then the
   runtime parameter :cpp:`fabarray.mfiter_tile_size` can be used to set it
   globally.

-  :cpp:`.validBox()` has the same meaning as in the non-tile approach,
   so we don’t use it. Instead, we use :cpp:`.tilebox()` to get the
   Box (and corresponding lo and hi) for the *current tile*, not the entire data region.

Let us consider an example. Suppose there are four boxes:

.. figure:: ./Parallel/domain-tile.png
   :alt: A simple domain showing 4 Boxes labeled 0–3, and their tiling regions (dotted lines)

   Boxes labeled 0–3, and their tiling regions (dotted lines)

The first box is divided into 4 logical tiles, the second and third
are divided into 2 tiles each (because they are small), and the fourth
into 4 tiles. So there are 12 tiles in total. The difference between
the tiling and non-tiling version are then:

-  In the tiling version,
   the loop body will be run 12 times. Note that :cpp:`tilebox` is
   different for each tile, whereas :cpp:`fab` might be referencing the
   same object if the tiles belong to the same box.

-  In the non-tiling
   version (by constructing MFIter without the optional second
   argument or setting to false), the loop body will be run 4 times
   because there are four boxes, and a call to :cpp:`mfi.tilebox()` will
   return the traditional validbox. The non-tiling case is
   essentially having one tile per box.

Tiling provides us the opportunity of a coarse-grained approach for
OpenMP. Threading can be turned on by inserting the following line
above the for (MFIter...) line.

::

      #pragma omp parallel

Assuming four threads are used in the above example, thread 0 will
work on 3 tiles from the first box, thread 1 on 1 tile from the first
box and 2 tiles from the second box, and so forth. Note that
OpenMP can be used even when tiling is turned off. In that case, the
OpenMP granularity is at the box level (and good performance would need
many boxes per MPI task).

While it is possible that, independent of whether or not tiling is on, OpenMP
threading could also be started within the function/loop called inside the
MFIter loop, rather than at the MFIter loop level, this
is not the approach taken in IAMR.

The tile size for the three spatial dimensions can be set by a
parameter, e.g., :cpp:`fabarray.mfiter_tile_size = 1024000 8 8`. A
huge number like 1024000 will turn off tiling in that direction.
As noted above, the MFIter constructor can also take an explicit
tile size: :cpp:`MFIter(mfi(mf,IntVect(128,16,32)))`.

Note that tiling can naturally transition from all threads working
on a single box to each thread working on a separate box as the boxes
coarsen (e.g., in multigrid).

The MFIter class provides some other useful functions:

::

     mfi.validbox()       : The same meaning as before independent of tiling.
     mfi.growntilebox(int): A grown tile box that includes ghost cells at box
                            boundaries only.  Thus the returned boxes for a
                            Fab are non-overlapping.
     mfi.nodaltilebox(int): Returns non-overlapping edge-type boxes for tiles.
                            The argument is for direction.
     mfi.fabbox()         : Same as mf[mfi].box().

Finally we note that tiling is not always desired or better. This
traditional fine-grained approach coupled with dynamic scheduling is
more appropriate for work with unbalanced loads, such as chemistry
burning in cells by an implicit solver. Tiling can also create extra
work in the ghost cells of tiles.

Practical Details in Working with Tiling
----------------------------------------

It is the responsibility of the coder to make sure that the routines within
a tiled region are safe to use with OpenMP. In particular, note that:

-  tile boxes are non-overlapping

-  the union of tile boxes completely cover the valid region of the fab

-  Consider working with a node-centered MultiFab, :cpp:`s_nd`, and a
   cell-centered MultiFab, :cpp:`s_cc`:

   -  with :cpp:`mfi(s_cc)`, the tiles are based on the cell-centered
      index space. If you have an :math:`8\times 8` box, then and 4 tiles, then
      your tiling boxes will range from :math:`0\rightarrow 3`, :math:`4\rightarrow
      7`.

   -  with :cpp:`mfi(s_nd)`, the tiles are based on nodal indices,
      so your tiling boxes will range from :math:`0\rightarrow 3`, :math:`4\rightarrow 8`.

-  When updating routines to work with tiling, we need to understand
   the distinction between the index-space of the entire box (which
   corresponds to the memory layout) and the index-space of the tile.
   Inside the MFIter loop, make sure to only update over
   the tile region, not for the entire box.
