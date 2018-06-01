# IAMR Regression Tests

Regression tests for IAMR are managed using the regression testing
tools that are distributed as part of the AMReX Codes, at the GitHub
repository "https://github.com/AMReX-Codes/regression_testing".  The
tools are a custom set of python scripts that check the results produced
by the code in the remote repository (i.e. not in your local, working
copy).  For information on how to set up a remote repository, see the
Setting Up the Tests section below.

The regression test scripts
1. "git pull" AMReX and the application code built on top of it, e.g.
   IAMR 
2. For each defined test,
   A. build an executable according to compile-time parameters/defines
   B. run in serial or parallel
   C. compare the results (typically in the form of a plotfile) with a
      "benchmark" reference solution.
3. Assemnble the results of all the tests, formatted in html.
   Optionally, if any of the tests fail, an email is generated and sent
   to a specified list of recipients.

A key design feature of the regression suite is that the reference
solutions can be updated manually at any time.  This is necessary
when, for example, a bug is discovered or an algorithm change results
in improved solutions or modified error metrics.  For this reason, the
initial set of benchmarks needs to be created manually before the first
test is executed (more info on how to do this in the Setting Up the
Tests section).


## Setting Up the Tests

1. A remote repository needs to be set up.  If that's not been done
already, go to IAMR on GitHub and fork the repo.  To do this, just
click on the "Fork" button located near the upper right hand corner of
the webpage.

2. Choose/create an area/directory to run the regression tests.  It is
recommended to make a special "scratch" area for the exclusive use of
the regression tester.  This is because the testing suite may optionally
checkout specific branches or SHA1 commits of the needed repositories,
and while it always restores the repository to it's original state, it's
best to let the tester work on its own where it won't risk confusing
you.

3. Clone the required repos into the scratch area.  Assuming the
following environment variable is set 

   *  REGTEST_SCRATCH: Scratch location set up in step 2 above

The following example commands will clone the required repositories
    ...
    git clone https://github.com/AMReX-Codes/amrex.git ${REGTEST_SCRATCH}/amrex
    git clone https://github.com/AMReX-Codes/regression_testing.git ${REGTEST_SCRATCH}/regression_testing
    git clone https://github.com/AMReX-Codes/IAMR.git ${REGTEST_SCRATCH}/IAMR
    ...

4.  Move to the location where the tests will be built/run, create
testing area expected by reg test scripts:

    ```
    cd ${REGTEST_SCRATCH}; mkdir -p TestData/IAMR
    ```

3.  Edit the config file, ${REGTEST_SCRATCH}/IAMR/Test/IAMR-tests.ini.
Set the AMReX and IAMR scratch clone locations and desired branch/SHA.
Also set

     testTopDir =  ${REGTEST_SCRATCH}/TestData/IAMR

4. Create symbolic links 
    ...
    ln -s ${REGTEST_SCRATCH}/regression_testing/regtest.py .
    ln -s ${REGTEST_SCRATCH}/IAMR/Test/IAMR-tests.ini .
    ...

4.  Generate the initial benchmark solutions for all the tests listed
in the .ini configuration file.  Rerunning this at any time will
overwrite the previous versions of the benchmarks

    ```
    ./regtest.py --make_benchmarks "<a useful comment>" IAMR-tests.ini
    ```

5. Upon some trigger event, re-run the tests and format the results in
html.  In this case, the results will appear as
TestData/IAMR/web/index.html

    ```
    ./regtest.py IAMR-tests.ini
    ```

NOTE: The regtest.py script takes a handy option "--no_update All",
which instructs the tester to work with the scratch repositories as
they currently exist.  Without this option specified, the branch of
each repository that is specified in the .ini config file is "git
pull"'d to obtain its most recent version; the original state of the
repositories are restored when the tests complete.  Using this
feature, a user can checkout any specific branch of any of the
repositories in the scratch area and run the complete set of tests.  A
user may wish to do this prior to issuing a "pull request", for
example.

More information on available options is given by
    ...
    ./regtest.py -h
    ...
which prints a verbose description of usage and setup. 