# IAMR Regression Tests

Regression tests for IAMR are managed using the regression testing
tools that are distributed as part of the AMReX Codes, at the GitHub
repository "AMReX-Codes/regression_testing".  The tools are a custom
set of python scripts that orchestrate repository "pull" operations of
AMReX and the application codes built on top of it.  For each defined
test, an executable is built according to compile-time
parameters/defines, and is run in serial or parallel; the results of
the runs (typically in the form of plotfiles) are compared with
"benchmark" reference solutions.  The results of the tests are
assembled, formatted in html.  Optionally, if any of the tests fail,
an email is generated and sent to a specified list of recipients.  A
key design feature of the regression suite is that the reference
solutions can be updated manually at any time.  This is necessary
when, for example, a bug is discovered or an algorithm change results
in improved solutions or modified error metrics.  For this reason, the
initial set of benchmarks need to be created manually before the first
test is executed.

This folder contains a configuration file for use with the AMReX
regression test suite that have been customized for IAMR. In order to
run the tests (either manually/locally or on an automated runner), a
scratch area on the runner machine is designated for each cloned
repository needed to build the tests, and the regression test results
(including the benchmark results).  NOTE: These locations are coded
into the test configuration file in IAMR's "Test" folder.

## Setting up the tests

The following example commands will clone the required repositories
into a scratch area, build the benchmarks and run the tests.  The
tests require a clone of the AMReX regression test suite, the AMReX
library, and the IAMR repository.  We assume that the following
environment variables are set for these steps

*  IAMR_REGTEST_HOME: Location where test code will be clone and run
*  AMREX_REGTEST_HOME: Location containing AMReX regression testing (cloned from 
`github.com:AMReX-Codes/regression_testing`)
*  AMREX_HOME: Location of scratch repository for AMReX
*  IAMR_HOME: Location of scratch `IAMR` repository to be tested


1.  Move to the location where the tests will be built/run, create
testing area expected by reg test scripts:

    ```
    cd ${IAMR_REGTEST_HOME}; mkdir -p TestData/IAMR
    ```

2.  Create a clone of `AMReX`, the AMReX `regression_testing` suite,
and `IAMR`.  Note that the testing suite may optionally checkout
specific branches or SHA1 commits of the needed repositories, but will
always restore the repository to it's original state afterward.
Nevertheless, it is suggested to make scratch clones of these
repositories for the exclusive use of the regression tester.

    ```
    git clone https://github.com/AMReX-Codes/amrex.git ${AMREX_HOME}
    git clone https://github.com/AMReX-Codes/regression_testing.git ${AMREX_REGTEST_HOME}
    git clone git@github.com:AMReX-Codes/IAMR.git ${IAMR_HOME}
    ```

3.  Generate the initial benchmark solutions for all the tests listed
in the .ini configuration file.  Rerunning this at any time will
overwrite the previous versions of the benchmarks

    ```
    ${AMREX_REGTEST_HOME}/regtest.py --make_benchmarks "<a useful comment>" ${IAMR_HOME}/Test/IAMR-tests.ini
    ```

4. Upon some trigger event, re-run the tests and format the results in
html.  In this case, the results will appear as
TestData/IAMR/www/index.html

    ```
    ${AMREX_REGTEST_HOME}/regtest.py ${IAMR_HOME}/Test/IAMR-tests.ini
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
