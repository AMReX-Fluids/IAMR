.. role:: cpp(code)
   :language: c++


.. _Chap:SetupAndRunning:

Creating and Running Your Own Problem
=====================================

This section covers how to set up your own problem and
provides some of the many options for controlling the algorithm and simulation runs.

In IAMR, a combination of coded routines and inputs parameters provide the problem
setup. Inputs parameters can also alter the behavior of the algorithm and
provide controls over simulation runs. IAMR uses AMReX's ``ParmParse`` class
infrastructure to read in the inputs parameters.
(See AMReX's documentation on ``ParmParse`` here: :ref:`sec:basics:parmparse`.)

Typically, the inputs parameters are collected into an inputs file, which is specified on the commandline.
Options can also be specified directly on the command line.

.. toctree::
   :caption: Contents:

   ProblemSetup
   RunningProblems
   AlgorithmOptions
