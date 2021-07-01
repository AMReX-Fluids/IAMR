.. role:: cpp(code)
   :language: c++

Debugging
===========================

Here we offer some tips we have found useful when debugging.

If the code is crashing, compile with ``DEBUG=TRUE`` and then examine the Backtrace file
produced from the running with the debug executable. Using the debug executable will also
turn on additional checks within the code (these are ``AMREX_ASSERT()`` statements).
Many times, these additional checks or the Backtrace
can point you to the problem. For more information and tips on debugging, also
see AMReX's documentation: :ref:`amrex:sec:basics:debugging`.

..
  If the issues is ``MLMG failed to converge`` 
