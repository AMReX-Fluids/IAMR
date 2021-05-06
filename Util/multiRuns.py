#!/usr/bin/env python3

# A template script to lunch several times IAMR at different resolution
# in order to evaluate the convergence order
# Used in the marex regression_testing framework

# Usage:
#   ./multiRuns.py --test_name DummyTest --input_file PeleLMInputFile

# Input:
#   * test_name: a TESTNAME that will be prependded to the plt files names
#   * PeleLMInputFile: the PeleLM input file

# "Internal" user input
#   * resolution : a list of the resolutions to run

# Head's up : 
#   * The PeleLM executable is searched for in the current directory.

import sys
import os
import shutil
import argparse
import numpy as np

USAGE = """    
    A template script to launch several times PeleLM.
"""

def multiRun(args):
    
    print(" Scripted multiple runs for convergence study ")
    # User data
    resolution = [32,64,128,256]

    # Get the PeleLM exec
    run_dir = os.getcwd()
    for f in os.listdir(run_dir):
        if ( f.startswith("amr2d") and f.endswith(".ex")):
               executable = f
    
    # Check the test name: current folder name is default
    if ( args.test_name == "None" ):
        args.test_name = run_dir.split("/")[-1]

    # Check for the input file: first inputs.* is default
    if ( args.input_file == "None" ):
        for f in os.listdir(run_dir):
            if ( f.startswith("inputs") ):
                args.input_file = f
                break

    # Loop on /= resolutions, run 
    for case in resolution:
        outfile = "{}_{}.run.out".format(args.test_name,case)
        print(" Running {}x{} case".format(case,case))
        runtime_params = "amr.n_cell={} {} {} ".format(case,case,case)
        runtime_params += "amr.plot_file={}_plt_{}_".format(args.test_name,case)
        os.system("mpiexec -n 1 ./{} {} {} > {}".format(executable, args.input_file, runtime_params, outfile))

def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument("--test_name", type=str, default="None", metavar="test-name",
                        help="name of the test. Default = current folder name")

    parser.add_argument("--input_file", type=str, default="None", metavar="Pele-input",
                        help="input file name. Default = first inputs.* in current directory")

    if not arg_string is None:
        args, unknown = parser.parse_known_args(arg_string)
    else:
        args, unknown = parser.parse_known_args()

    return args   

if __name__ == "__main__":
    args = parse_args(arg_string=sys.argv[1:])
    multiRun(args)
