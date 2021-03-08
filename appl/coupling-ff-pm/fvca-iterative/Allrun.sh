#! /usr/bin/env bash

solver_input="params.input"

ff_solver="test_ff_reversed"
pm_solver="test_pm_reversed"

precice_config="precice-config-serial-implicit-reversed.xml"

rm -rf "precice-run/"
ff_cmd="./${ff_solver} - ${precice_config}"
echo "${ff_cmd}"
./${ff_solver} - ${precice_config} > ${ff_solver}.log 2>&1 &
PIDFluid=$!
./${pm_solver} - ${precice_config} > ${pm_solver}.log 2>&1 &
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
