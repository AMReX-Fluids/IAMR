AMR
=====

Our approach to adaptive refinement in IAMR uses a nested
hierarchy of logically-rectangular grids with simultaneous refinement
of the grids in both space and time. The integration algorithm on the grid hierarchy
is a recursive procedure in which coarse grids are advanced in time,
fine grids are advanced multiple steps to reach the same time
as the coarse grids and the data at different levels are then synchronized.

During the regridding step, increasingly finer grids
are recursively embedded in coarse grids until the solution is
sufficiently resolved. An error estimation procedure based on
user-specified criteria (described in :ref:`Tagging for Refinement`)
evaluates where additional refinement is needed
and grid generation procedures dynamically create or
remove rectangular fine grid patches as resolution requirements change.

The dynamic creation and destruction of grid levels is a fundamental part of IAMR’s capabilities.
At regular intervals (set by the user), each Amr level that is not the finest allowed for the run
will invoke a “regrid” operation. When invoked, a set of error tagging functions is traversed. For each,
a field specific to that function is derived from the state over the level, and passed through a kernel
that “set”’s or “clear”’s a flag on each cell.
The field and function for each error tagging quantity is
identified in the setup phase of the code (in NS\_error.cpp).
Each function adds or removes to the list of cells tagged for refinement.
This may then be extended if amr.n\_error\_buf :math:`> 0` to a certain number
of cells beyond these tagged cells.
This final list of tagged
cells is sent to a grid generation routine, which uses the Berger-Rigoutsos algorithm to create rectangular grids
which will define a new finer level (or set of levels). State data is filled over these new grids, copying where
possible, and interpolating from coarser level when no fine data is available. Once this process is complete,
the existing Amr level(s) is removed, the new one is inserted into the hierarchy, and the time integration
continues.

Tagging for Refinement
---------------------------------------

This section describes the process
by which cells are tagged during regridding, and describes how to add customized tagging criteria.
IAMR provides two methods for creating error estimation functions: dynamic generation and explicit or
hard-coded. First dynamic generation is discussed, followed by how to explicitly define error functions.

Dynamically generated tagging functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dynamically created error functions are based on runtime data specified in the
inputs (ParmParse) data.
These dynamically generated error functions can tag a specified region of the domain and/or test on a field
that is a state variable
or derived variable defined in NS\_derive.cpp and included in the derive\_lst in NS\_setup.cpp.
Available tests include

-  “greater\_than”: :math:`field >= threshold`

-  “less\_than”: :math:`field <= threshold`

-  “adjacent\_difference\_greater”: :math:`max( | \text{difference between any nearest-neighbor cell} | ) >= threshold`

-  “vorticity”: :math:` |vorticity| >= 2^{level} * threshold`

The following example portions of ParmParse’d input files demonstrate the usage of this feature.
This first example tags all cells inside the region ((.25,.25)(.75,.75)):

::

          amr.refinement_indicators = box

          amr.box.in_box_lo = .25 .25
          amr.box.in_box_hi = .75 .75

The next example adds three user-named criteria –
hi\_rho: cells with density greater than 1 on level 0, and greater than 2 on
level 1 and higher;
lo\_temp: cells with T less than 450K that are inside the region ((.25,.25)(.75,.75));
and Tdiff: cells having a temperature difference of 20K
or more from that of their
immediate neighbor. The first will trigger up to Amr level 3, the second only to level 1, and the third to level 2.
The third will be active only when the problem time is between 0.001 and 0.002 seconds.

::

          amr.refinement_indicators = hi_rho lo_temp Tdiff

          amr.high_rho.max_level = 3
          amr.hi_rho.value_greater = 1. 2.
          amr.hi_rho.field_name = density

          amr.lo_temp.max_level = 1
          amr.lo_temp.value_less = 450
          amr.lo_temp.field_name = temp
          amr.lo_temp.in_box_lo = .25 .25
          amr.lo_temp.in_box_hi = .75 .75

          amr.Tdiff.max_level = 2
          amr.Tdiff.adjacent_difference_greater = 20
          amr.Tdiff.field_name = temp
          amr.Tdiff.start_time = 0.001
          amr.Tdiff.end_name = 0.002

Note that these criteria can be modified between restarts of the code.
By default, the new criteria will take effect at the next
scheduled regrid operation. Alternatively, the user may restart with amr.regrid\_on\_restart = 1 in order to
do a full (all-levels) regrid after reading the checkpoint data and before advancing any cells.

Explicitly defined tagging functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Explicitly defined error estimation functions can be used either instead of or in addition to
dynmaically generated funtions. These functions can be added to NavierStokes::errorEst() in
NS\_error.cpp. Any dynamically generated error functions will operate first.
Please note that while CLEARing a tagged cell is possible, it is not reccomended as it
may not have the desired effect.
