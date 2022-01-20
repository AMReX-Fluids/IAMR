.. role:: cpp(code)
   :language: c++

Debugging
===========================

Here we offer some tips we have found useful when debugging.

If the code is crashing, compile with ``DEBUG=TRUE`` and then examine the Backtrace file
produced from the running with the debug executable. Using the debug executable will also
turn on additional checks within the code (these are ``AMREX_ASSERT()`` statements).
Many times, these additional checks or the Backtrace can point you to the problem.

It can be helpful to simplify the problem as much as possible before digging into debugging.
For example, reduce the problem size and set ``amr.max_level`` as small as possible.
Diffusion and viscosity and be turned off by setting ``ns.vel_visc_coef = 0`` and
``ns.scal_diff_coefs = 0``.

Sometimes, restarting from a checkpoint file can reduce the amount of time spent waiting for
the code to fail durning the debugging process.

If the problem only happens with multilevel runs, ``amr.subcycling_mode = None`` will turn
time subcycling off and advance the whole system at the same dt. 

Please also read AMReX's documentation: :ref:`amrex:sec:basics:debugging`
for more information and tips.

If you believe you've encountered a bug or incorrect behavior in IAMR, please report the issue
on IAMR's github page `here <https://github.com/AMReX-Codes/IAMR/issues>`_ .

..
  If the issues is ``MLMG failed to converge`` 
