
Controlling What’s in the PlotFile
==================================

There are a few options that can be set at runtime to control what
variables appear in the plotfile. See section :ref:`sec:InputsPlotfiles`.

Amrvis
======

Amrvis is a package developed
by CCSE that is designed specifically for highly efficient visualization
of block-structured hierarchical AMR data.
Information on building the amrvis2d and amrvis3d executables is given in section
:ref:`Visualizing the Results`. A very useful feature of Amrvis is
View :math:`\rightarrow` Dataset, which
allows you to actually view the numbers in a spreadsheet that is nested
to reflect the AMR hierarchy – this can be handy for
debugging. You can modify how many levels of data you want to see,
whether you want to see the grid boxes or not, what palette you use,
etc.

If you like to have amrvis display a certain variable, at a certain
scale, when you first bring up each plotfile (you can always change it
once the amrvis window is open), you can modify the amrvis.defaults
file in your directory to have amrvis default to these settings every
time you run it. The Amrvis repository has a sample amrvis.defaults
file that can be copied to your local folder and modified.

VisIt
=====

VisIt reads AMReX plotfiles, including IAMR's, natively.  Open the ``Header`` file inside
the plotfile directory (e.g. ``plt00000/Header``) with "File" :math:`\rightarrow` "Open file";
VisIt picks its Boxlib reader for it automatically.  See the :ref:`VisIt section of AMReX's
documentation <amrex:sec:visit>` for a worked example, and
https://visit-dav.github.io/visit-website/ for VisIt itself.

yt
=====

yt (https://yt-project.org) is an open-source Python package for analyzing and visualizing
simulation data, and it reads IAMR plotfiles directly; IAMR's fields appear under the
``"boxlib"`` field type, e.g. ``("boxlib", "density")``.  Install it with
``python -m pip install yt``.  For how to use it with AMReX plotfiles, see the
:ref:`yt section of AMReX's documentation <amrex:section-2>` and the
`AMReX section <https://yt-project.org/doc/examining/loading_data.html#amrex-boxlib-data>`_
of yt's guide to loading data.
