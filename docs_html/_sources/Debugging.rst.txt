.. role:: cpp(code)
   :language: c++

.. _Chap:Debugging:

Debugging
=========

Debugging is an art.  Everyone has their own favorite method.  Here we
offer a few tips we have found to be useful.

Compiling in debug mode (e.g., :cpp:`make DEBUG=TRUE` for gmake users;
:cpp:`cmake -DDEBUG` for cmake users) and running with
:cpp:`amrex.fpe_trap_invalid=1` in the inputs file can be helpful.
In debug mode, many compiler debugging flags are turned on and all
:cpp:`MultiFab` s are initialized to signaling NaNs.  The
:cpp:`amrex.fpe_trap_invalid` parameter will result in backtrace files
when a floating point exception occurs.  One can then examine those
files to track down the origin of the issue.

Several other ways to look at the data include:

1) Writing a :cpp:`MultiFab` to disk with

.. highlight:: c++

::

    VisMF::Write(const FabArray<FArrayBox>& mf, const std::string& name);

and examining it with ``Amrvis`` (section :ref:`sec:amrvis` in the AMReX documentation).

2) You can also use the :cpp:`print_state` routine: 

.. highlight:: c++

::

    void print_state(const MultiFab& mf, const IntVect& cell, const int n=-1);

which outputs the data for a single cell.

3) If you want to compare old and new plotfiles, 

.. highlight:: c++

::

    fcompare --infile1 plt00000_run1 --infile2 plt00000_run2 --diffvar u_g

will print out the maximum absolute and relative differences between the two plotfiles
for each variable and will also create a new plotfile "diffs" that contains the difference
in u_g (in this case) between the two plotfiles.

The :cpp:`fcompare` executable can be built in AMReX (go to amrex/Tools/Plotfile and type "make").

4) Valgrind is another useful debugging tool.  Note that for runs using
more than one MPI process, one can tell valgrind to output to different 
files for different processes.  For example,

.. highlight:: console

::

    mpiexec -n 4 valgrind --leak-check=yes --track-origins=yes --log-file=vallog.%p ./amr3d.exe ...
