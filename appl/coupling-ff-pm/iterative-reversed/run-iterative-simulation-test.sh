#! /usr/bin/env bash


ff_solver="test_ff_reversed"
pm_solver="test_pm_reversed"


rm -f *.csv *.pvd *.vtu *.log *.txt *.json

#echo "" > ${configurations}
basedir="${PWD}"

if [ $# -ne 2 ]; then
    echo "Number of parameters is wrong! (${#} instead of 3)"
    exit 1
fi


dumux_param_file=$1
precice_config_file=$2

echo "Dumux parameter file: ${dumux_param_file}"
echo "preCICE configuration file: ${precice_config_file}"
            
# Check if Stokes or Navier-Stokes        
flowProblemName="stokes"
#          echo ${hasInertiaTerms}
if [[ "${hasInertiaTerms}" == "true" ]]; then
    flowProblemName="navier-stokes"
fi


rm -rf "precice-run/"
ff_cmd="./${ff_solver} ${dumux_param_file} - ${precice_config_file}"
echo "${ff_cmd}"
./${ff_solver} ${dumux_param_file} - ${precice_config_file} > ${ff_solver}.log 2>&1 &
PIDFluid=$!
./${pm_solver} ${dumux_param_file} - ${precice_config_file} > ${pm_solver}.log 2>&1 &
PIDSolid=$!

echo "Waiting for the participants to exit..."
echo "(you may run 'tail -f ${ff_solver}.log' or 'tail -f ${pm_solver}.log' in another terminal to check the progress)"

wait ${PIDFluid}
wait ${PIDSolid}

if [ $? -ne 0 ] || [ "$(grep -c -E "error:" ${ff_solver}.log)" -ne 0 ] || [ "$(grep -c -E "error:" ${pm_solver}.log)" -ne 0 ]; then
    echo ""
    echo "Something went wrong... See the log files for more."
else
    echo ""
    echo "The simulation completed!"
fi






