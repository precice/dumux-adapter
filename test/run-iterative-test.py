#!/usr/bin/env python3

import argparse
import json
import subprocess as sp
import sys
import shutil

# Make sure that the PYTHONPATH is set properly
from fuzzycomparevtu import compare_vtk



# parse arguments
parser = argparse.ArgumentParser()
#parser.add_argument('-c', '--command', nargs=1, help='The executable and optional arguments as a single string', required=True)
parser.add_argument('-f', '--files', nargs='+', help="Pairs of file names (first reference, then current). Usage: '[-f ref1 cur1 [[ref2] [cur2] ...]]'")
parser.add_argument('-r', '--relative', type=float, default=1e-2, help='maximum relative error (default=1e-2) when using fuzzy comparison')
parser.add_argument('-a', '--absolute', type=float, default=1.5e-7, help='maximum absolute error (default=1.5e-7) when using fuzzy comparison')
parser.add_argument('-z', '--zeroThreshold', type=json.loads, default='{}', help='Thresholds for treating numbers as zero for a parameter as a python dict e.g. {"vel":1e-7,"delP":1.0}')
parser.add_argument('-dpf', '--dumux-param-file', type=str, default='params.input', help='DuMuX parameter file (required)', required=True)
parser.add_argument('-pcf', '--precice-config-file', type=str, default='precice-config.xml', help='preCICE configuration file (required)', required=True)
args = vars(parser.parse_args())
print(args)


def run_solver( solver_name, dumux_param_file, precice_config_file ):
    print( "Starting {}".format(solver_name) )
    f = open("{}.log".format(solver_name), 'w')

    #print( "./{}".format(solver_name) )
    #print( str(dumux_param_file[0]) )
    #print( precice_config_file )

    proc = sp.Popen(["./{}".format(solver_name), dumux_param_file, "-", precice_config_file], stdout=f, stderr=f)
    return proc, f

try:
    shutil.rmtree("precice-run/")
except:
    print("no directory \"precice-run/\" to delete")

# Start thread for Biot solver in background
ff_proc, ff_output = run_solver("test_ff_reversed", args['dumux_param_file'], args['precice_config_file'] )
pm_proc, pm_output = run_solver("test_pm_reversed", args['dumux_param_file'], args['precice_config_file'] )

# Wait for solvers to finish
pm_proc.wait()
ff_proc.wait()


return_code = 0
for i in range(0, len(args['files'])//2):
    print("\nFuzzy comparison...")
    result = compare_vtk(args['files'][i*2], args['files'][(i*2)+1], relative=args['relative'], absolute=args['absolute'], zeroValueThreshold=args['zeroThreshold'])
    if result:
        return_code = 1

if return_code==0:
    try:
        shutil.rmtree("precice-run/")
    except:
        print("no directory \"precice-run/\" to delete")

sys.exit(return_code)

