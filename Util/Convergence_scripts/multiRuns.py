#!/usr/bin/env python3

# A template script to launch IAMR several times at different resolutions
# in order to evaluate the convergence order

# Usage:
#   ./multiRuns.py --test_name DummyTest --exe_name Executable --input_file InputFile

# Input:
#   * test_name: a TESTNAME that will be prependded to the plt files names
#   * Executable: the IAMR executable
#   * InputFile: the input file

# "Internal" user input
#   * resolution : a list of the resolutions to run

# Head's up : 
#   * The executable is searched for in the current directory.

import sys
import os
import shutil
import argparse
import numpy as np

USAGE = """    
    A template script to launch IAMR several times
"""

def multiRun(args):
    
    print(" Scripted multiple runs for convergence study ")
    # User data
    resolution = [32,64,128,256]

    # Get the current directory
    run_dir = os.getcwd()
    
    # Get the executable: first amr*.ex is default
    if ( args.exe_name == "None" ):
        for f in os.listdir(run_dir):
            if ( f.startswith("amr") and f.endswith(".ex")):
                executable = f
                break
    else:
        executable = args.exe_name

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
        print(" Running {} on {}x{} case".format(executable,case,case))
        runtime_params = "amr.n_cell={} {} {} ".format(case,case,case)
        runtime_params += "amr.plot_file={}_plt_{}_".format(args.test_name,case)
        os.system("mpiexec -n 2 ./{} {} {} > {}".format(executable, args.input_file, runtime_params, outfile))

def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument("--test_name", type=str, default="None", metavar="test-name",
                        help="name of the test. Default = current folder name")

    parser.add_argument("--exe_name", type=str, default="None", metavar="exe-name",
                        help="name of the executable. Default = first amr*.ex in current directory")

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
