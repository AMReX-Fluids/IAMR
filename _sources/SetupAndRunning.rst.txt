.. role:: cpp(code)
   :language: c++


.. _Chap:SetupAndRunning:

Creating and Running Your Own Problem
=====================================

This section covers how to set up your own problem and
provides some of the many options for controlling the algorithm and simulation runs.

In IAMR, a combination of coded routines and inputs parameters provide the problem
setup. Inputs parameters can also alter the behavior of the algorithm and
provide controls over simulation runs.
Typically, the inputs parameters are collected into an inputs file, which is specified on
the commandline. Options can also be specified directly on the command line.
IAMR uses AMReX's ``ParmParse`` class
infrastructure to read in the inputs parameters.
(See AMReX's documentation on ``ParmParse`` here: :ref:`sec:basics:parmparse`.)

.. tip::
   Before making any code changes, we highly reccomend creating a new git branch.

Extra points if the new branch's name reflects the changes, e.g. ``couette_flow``.
We do `not` reccomend making any changes within the ``development`` branch because
doing so will make your local ``development`` branch diverge from the matching branch
in the main IAMR repository. This can lead to difficulties pulling in updates and bug-fixes
from the main IAMR repo.

Here, we draw your attention to the option of creating a github "fork" of IAMR.
A fork is a copy of a repository that you manage (see the github documentation for more details
on forks).
For more information on creating a fork of IAMR, see the :ref:`Chap:Contributing` Section (and yes,
this is a reasonable option even if you don't plan to contribute code changes to IAMR).

If you choose not to create a fork, you can create a new git branch in your local clone by
executing the following commands within the ``IAMR`` directory:

.. code:: shell

  git checkout development
  git checkout -b <branch_name>

To pull in any updates from the IAMR repo, execute:

.. code:: shell

  git fetch
  git merge development


.. toctree::
   :caption: Contents:

   ProblemSetup
   RunningProblems
   AlgorithmOptions
