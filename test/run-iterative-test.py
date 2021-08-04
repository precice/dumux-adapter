#!/usr/bin/env python3

import argparse
import glob
import json
import os
import subprocess as sp
import sys
import shutil

# Make sure that the PYTHONPATH is set properly
from fuzzycomparevtu import compare_vtk


# parse arguments
parser = argparse.ArgumentParser()
# parser.add_argument('-c', '--command', nargs=1, help='The executable and optional arguments as a single string', required=True)
parser.add_argument(
    "-f",
    "--files",
    nargs="+",
    help="Pairs of file names (first reference, then current). Usage: '[-f ref1 cur1 [[ref2] [cur2] ...]]'",
)
parser.add_argument(
    "-r",
    "--relative",
    type=float,
    default=1e-2,
    help="maximum relative error (default=1e-2) when using fuzzy comparison",
)
parser.add_argument(
    "-a",
    "--absolute",
    type=float,
    default=1.5e-7,
    help="maximum absolute error (default=1.5e-7) when using fuzzy comparison",
)
parser.add_argument(
    "-z",
    "--zeroThreshold",
    type=json.loads,
    default="{}",
    help='Thresholds for treating numbers as zero for a parameter as a python dict e.g. {"vel":1e-7,"delP":1.0}',
)
parser.add_argument(
    "-dpf",
    "--dumux-param-file",
    type=str,
    default="params.input",
    help="DuMuX parameter file (required)",
    required=True,
)
parser.add_argument(
    "-pcf",
    "--precice-config-file",
    type=str,
    default="precice-config.xml",
    help="preCICE configuration file (required)",
    required=True,
)
parser.add_argument(
    "-cn",
    "--case-name",
    type=str,
    #    default="precice-config.xml",
    help="Name identifier of the test case (required)",
    required=True,
)
parser.add_argument(
    "-pif",
    "--precice-iteration-files",
    type=str,
    #    default="precice-config.xml",
    nargs=2,
    help="Name of the preCICE iteration file names to be compared (required)",
    required=True,
)
args = vars(parser.parse_args())
print(args)


def run_solver(solver_name, dumux_param_file, precice_config_file):
    print("Starting {}".format(solver_name))
    f = open("{}.log".format(solver_name), "w")

    # print( "./{}".format(solver_name) )
    # print( str(dumux_param_file[0]) )
    # print( precice_config_file )

    proc = sp.Popen(
        ["./{}".format(solver_name), dumux_param_file, "-", precice_config_file],
        stdout=f,
        stderr=f,
    )
    return proc, f


def diff_iteration_files(diff_file_name, file_names):
    assert len(file_names) == 2, "Script expects two iteration files to be compared"

    print(
        "Start diff on files:\n  File 1: {}\n  File 2: {}".format(
            file_names[0], file_names[1]
        )
    )

    with open(diff_file_name, "w") as f:
        proc = sp.Popen(
            ["diff", file_names[0], file_names[1]],
            stdout=f,
            stderr=f,
        )
        proc.wait()
        return_code = proc.returncode

    if return_code == 0:
        print("Diff succeeded.")
    else:
        print("Diff failed. Return code: {}".format(return_code))

    return return_code


try:
    shutil.rmtree("precice-run/")
except:
    print('no directory "precice-run/" to delete')

# Start thread for Biot solver in background
ff_proc, ff_output = run_solver(
    "fvca-iterative-ff", args["dumux_param_file"], args["precice_config_file"]
)
pm_proc, pm_output = run_solver(
    "fvca-iterative-pm", args["dumux_param_file"], args["precice_config_file"]
)

# Wait for solvers to finish
pm_proc.wait()
ff_proc.wait()

print("Comparing VTU files")
return_code = 0
for i in range(0, len(args["files"]) // 2):
    print("\nFuzzy comparison...")

    print(
        "\n  Compare: \n    File 1: {}\n    File 2: {}".format(
            args["files"][i * 2], args["files"][i * 2 + 1]
        )
    )
    result = compare_vtk(
        args["files"][i * 2],
        args["files"][(i * 2) + 1],
        relative=args["relative"],
        absolute=args["absolute"],
        zeroValueThreshold=args["zeroThreshold"],
    )
    if result:
        return_code = 1

print("Comparing preCICE iteration files")
if (
    diff_iteration_files(
        "{case_name}-diff.txt".format(case_name=args["case_name"]),
        args["precice_iteration_files"],
    )
    != 0
):
    return_code = 1


print("Moving relevant files")
if return_code != 0:
    try:
        shutil.move(
            "fvca-iterative-ff.log",
            "{case_name}_ff.log".format(case_name=args["case_name"]),
        )
        shutil.move(
            "fvca-iterative-pm.log",
            "{case_name}_pm.log".format(case_name=args["case_name"]),
        )
    except:
        print("Not all files to be saved were found!")

print("Removing files that are not needed anymore")
if return_code == 0:
    try:
        shutil.rmtree("precice-run/")
        # shutil.rm("*.pvd")
        # shutil.rm("*.log")
        # shutil.rm("*.json")
    except:
        print('no directory "precice-run/" to delete')

    try:
        files_to_delete = glob.glob("./*.vtu")
        files_to_delete += glob.glob("./*.pvd")
        files_to_delete += glob.glob("./*.json")
        files_to_delete += glob.glob("./*.log")
        files_to_delete += glob.glob("./*.txt")
        print(files_to_delete)

        # os.remove( files_to_delete )
        for file_to_delete in files_to_delete:
            os.remove(file_to_delete)
        #    shutil.rm( file_to_delete )
        # shutil.rm( glob.glob("./*.pvd"))
        # shutil.rm( glob.glob("./*.json"))
    except:
        print("No further files to remove")

sys.exit(return_code)
