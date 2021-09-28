#!/usr/bin/env python3

# Template post-processing script for IAMR convergence analysis
# Must be used after multirun.py script
# Input are limited by the regression framework.

# Usage:
#   ./pprocConvOrder.py --pproc_exe prog.exe --test_name DummyTest

# Input:
#   * --pproc_exe: the processing executable path
#   * --test_name: a TESTNAME that will be looked for during postprocessing


# TODO -- add these command arguments. Use python to compile diff tool internally.
# combine into one script with multiruns.py.
# Usage:
#   ./pprocConvOrder.py --diff_tool tool_name --vars "var1 var2" --resolution "res1 res2 res3 res4" --test_name DummyTest

# Input:
#   * --diff_tool : the AMReX plotfile differencing tool. Can be either fcompare or diffsamedomain.
#                   fcompare: difference plotfiles with the same resolution. Here it's used to get the error from the initial solution (== analytical solution)
#                   diffsamedomain: difference plotfiles with the same domain but different resolutions. Analytical solution is not known and errors are computed from the next finer grid
#   * --vars : a list of the variables of interest (no check is done on whether it exists in the plotfile ...)
#   * --resolution : a list of the resolutions to post-process (should be consistent with multirun.py, if used)
#   * --test_name : a TESTNAME that will be looked for during postprocessing


# Output:
#   * Convergence_${TESTNAME}.png file with the log-log plot of the error vs. resolution.
#   * ConvTable_${TESTNAME}.tex file with the convergence rate formatted in an LaTeX table.
#   * Convergence_${TESTNAME}.dat plain text file with the convergence rate.

# Head's up : 
#   - The script will get a copy of the post-processing program (if not already there) in the testing folder. The name of this folder is assumed to be the TESTNAME.  
#   - The plt files naming convention is: ${TESTNAME}_plt_${resolution}_*****. It is used to get the first and last solution of a test at a given resolution.
#   - Errors are parsed from the screen output of the standard fcompare/diffsamedomain. Beware of any change of these programs. 

import sys
import os
import fnmatch
import shutil
import argparse
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

USAGE = """
    Template post-processing script for convergence analysis
"""

def pproc(args):

    # User data
    vars=["y_velocity","x_velocity"]
    resolution = [32,64,128,256]
    pproc_type = "fcompare"

    # Get a local copy of post-processing executable
    run_dir = os.getcwd()
    if ( not os.path.isfile(os.path.basename(args.pproc_exe)) ):
        shutil.copy(args.pproc_exe, run_dir)

    # Check the test name: current folder name is default
    if ( args.test_name == "None" ):
        args.test_name = run_dir.split("/")[-1]

    # Run the postprocessing
    if ( pproc_type == "fcompare" ):     # running fcompare since analytical solution is known
        errors = np.empty([len(resolution),len(vars)+1])
        pltfile=[]
        for res in range(len(resolution)):
            case = resolution[res]
            errors[res,0] = case

            # Get the fcompare inputs: first and last solution of current case
            # TODO: the analytical solution might not be plt****_00000 ...
            for f in os.listdir(run_dir):
                if ( not fnmatch.fnmatch(f, '*old*')):
                    if (f.startswith("{}_plt_{}_".format(args.test_name,case))):
                        pltfile.append(f)
            pltfile.sort()
            outfile = "error_{}.analysis.out".format(case)

            print("running: ./{} -a -n 2 {} {}".format(os.path.basename(args.pproc_exe), pltfile[0], pltfile[-1])) 

            os.system("./{} -a -n 2 {} {} > {}".format(os.path.basename(args.pproc_exe), pltfile[0], pltfile[-1], outfile))
            pltfile.clear()
        
            # Extract errors on each variable
            with open(outfile) as fp:
                for i, line in enumerate(fp):
                    if (i >= 5): #throw away first 5 lines, these are header
                        var = line.split()[0]
                        for v in range(len(vars)):
                            if ( var == vars[v] ):
                                errors[res,v+1] = line.split()[1]
            os.system("rm {}".format(outfile))
    elif ( pproc_type == "diffsamedomain" ):   # running diffsamedomain. No analytical sol ...
        errors = np.empty([len(resolution)-1,len(vars)+1])
        pltfile=[]
        pltfilenext=[]
        for res in range(len(resolution)-1):
            case = resolution[res]
            nextcase = resolution[res+1]
            errors[res,0] = case

            # Get the diffsamedomain inputs: last solutions of current 
            # and next finer cases. These run should have been runned to the same final time
            for f in os.listdir(run_dir):
                if ( not fnmatch.fnmatch(f, '*old*')):
                    if (f.startswith("{}_plt_{}_".format(args.test_name,case))):
                        pltfile.append(f)
                    if (f.startswith("{}_plt_{}_".format(args.test_name,nextcase))):
                        pltfilenext.append(f)
            pltfile.sort()
            pltfilenext.sort()

            print("command: ./{} -a -n 2 {} {}".format(os.path.basename(args.pproc_exe), pltfile[-1], pltfilenext[-1])) 

            outfile = "error_{}.analysis.out".format(case)
            os.system("./{} infile1={} reffile={} > {}".format(os.path.basename(args.pproc_exe), pltfile[-1], pltfilenext[-1], outfile))
            pltfile.clear()
            pltfilenext.clear()

            # Extract errors on each variable
            with open(outfile) as fp:
                for i, line in enumerate(fp):
                    if (i >= 4):
                        var = line.split(":")[0]
                        for v in range(len(vars)):
                            if ( var.split(" ")[0] == vars[v] ):
                                errors[res,v+1] = line.split(":")[1]
            os.system("rm {}".format(outfile))
    else:
        print("Wrong pproc_type: {}. should be either fcompare or diffsamedomain".format(pproc_type))
        return

    print("\nRaw error data: \n{}".format(errors))

    leastSq = np.empty(len(vars))
    
    # Plot data
    plotdata(errors, args.test_name, vars, leastSq)
    # Write out data in a tex file 
    writetex(errors, args.test_name, vars, leastSq)
    # Write out data as a txt file
    writeRegTestFile(errors, args.test_name, vars, leastSq)

def plotdata(data, test_name, vars, leastSq):
    # Compute convergence order with least squares fit
    x = data[:,0]
    logx=np.log10(x)
    print("Convergence order computed via least squares fit:")
    for i in range(0, len(vars)):
        y = data[:,i+1]
        logy = np.log10(y)
        #plt.plot(x, y, 'o', label='Original data', markersize=10)
        A = np.vstack([logx, np.ones(len(x))]).T
        m, c = np.linalg.lstsq(A, logy, rcond=None)[0]
        print("  {} {}".format(vars[i], m))
        plt.plot(x, x**m * 10**c, linestyle='dotted', label="Fitted line for {}".format(vars[i]))
        leastSq[i] = m
        
    # Evaluate 2nd order slope
    snd_order = data[:,1]*1.05
    for i in range(1,len(data[:,1])):
        snd_order[i] = snd_order[i-1]/np.exp(2.0*np.log(2.0))
    for i in range(0, len(vars)):
        plt.plot(data[:,0], data[:,i+1], label="{}".format(vars[i]))
    plt.plot(data[:,0], snd_order[:],linestyle='--',color='k', label='2nd-order')
    plt.xlabel("Resolution")
    plt.ylabel("Error L2norm")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(which='both',color='k', linestyle=':', linewidth=1)
    plt.legend(bbox_to_anchor=(0.9, 0.9), loc=1, borderaxespad=0.)
    plt.savefig("Convergence_{}.png".format(test_name))

def writetex(data, test_name, vars, leastSq):
    # Evaluate order
    conv_order = np.empty([len(data[:,0])-1,len(vars)])
    for v in range(len(vars)):
        for i in range(len(conv_order[:,0])):
            conv_order[i,v] = np.log(data[i,v+1]/data[i+1,v+1])/np.log(2.0)
    fout = open("ConvTable_{}.tex".format(test_name), "w")            
    fout.write("\\begin{table}[ht!]\n")
    fout.write("\centering\n")
    fout.write("\\begin{tabular}{l|")
    for i in range(len(conv_order[:,0])):
        fout.write("c ")
    fout.write("}\n")    
    fout.write("\hline\n")
    fout.write("Variable ")
    for i in range(len(conv_order[:,0])):
        fout.write("&  {}/{} ".format(data[i+1,0],data[i,0]))
    fout.write("&  least squares fit ")
    fout.write("\\\\\n\hline\hline\n")
    for v in range(len(vars)):
        fout.write("{} ".format(vars[v].replace("_","\_")))
        for i in range(len(conv_order[:,0])):
            fout.write("&  {:.3f} ".format(conv_order[i,v]))
        fout.write("&  {:.3f} ".format(leastSq[v]))
        fout.write("\\\\\n")
    fout.write("\end{tabular}\n")
    fout.write("\caption{convergence order}\n")
    fout.write("\label{table:conv}\n")
    fout.write("\end{table}\n")
    fout.close()

def writeRegTestFile(data, test_name, vars, leastSq):
    # Evaluate order
    conv_order = np.empty([len(data[:,0])-1,len(vars)])
    for v in range(len(vars)):
        for i in range(len(conv_order[:,0])):
            conv_order[i,v] = np.log(data[i,v+1]/data[i+1,v+1])/np.log(2.0)
    fout = open("Convergence_{}.dat".format(test_name), "w")            
    fout.write(" Variables ")
    for i in range(len(conv_order[:,0])):
        fout.write(" {}/{} ".format(data[i+1,0],data[i,0]))
    fout.write(" least squares fit \n")
    for v in range(len(vars)):
        fout.write("{} ".format(vars[v]))
        for i in range(len(conv_order[:,0])):
            fout.write("  {:.3f} ".format(conv_order[i,v]))
        fout.write("  {:.3f} ".format(leastSq[v]))
        fout.write("\n")
    fout.close()

def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument("--test_name", type=str, default="None", metavar="test-name",
                        help="name of the test. Default = current folder name")

    parser.add_argument("--pproc_exe", type=str, default="None", metavar="pproc.exe",
                        required=True,help="path to the executable required for the analysis.")

    if not arg_string is None:
        args, unknown = parser.parse_known_args(arg_string)
    else:
        args, unknown = parser.parse_known_args()

    return args   

if __name__ == "__main__":
    args = parse_args(arg_string=sys.argv[1:])
    pproc(args)
