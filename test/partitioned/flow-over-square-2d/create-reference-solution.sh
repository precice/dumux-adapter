#! /usr/bin/env bash

# Start progam
#solver_input="params.input"

ff_solver="fvca-iterative-ff"
pm_solver="fvca-iterative-pm"


function clean_up_dir()
{
    echo "Cleaning up directory"
    rm -rf "precice-run/"
    rm -rf *.json
    rm -rf *.pvd
    rm -rf *.vtu
}


function run_test_case()
{
    #precice_config="precice-config-serial-implicit-reversed.xml"
    if [[ $# -eq 2 ]]; then
        solver_input="${1}"
        precice_config="${2}"
    else
        echo "ERROR: No preCICE configuration was specified."
        exit 1
    fi

    echo "Running case using ${solver_input} - ${precice_config}"

    clean_up_dir
    ff_cmd="./${ff_solver} ${solver_input} - ${precice_config}"
    echo "${ff_cmd}"
    #./${ff_solver} - ${precice_config} > ${ff_solver}.log 2>&1 &
    ${ff_cmd} > ${ff_solver}.log 2>&1 &
    PIDFluid=$!
    pm_cmd="./${pm_solver} ${solver_input} - ${precice_config}"
    #./${pm_solver} - ${precice_config} > ${pm_solver}.log 2>&1 &
    ${pm_cmd} > ${pm_solver}.log 2>&1 &
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
}

function move_result_to_dir()
{
    if [[ $# -eq 2 ]]; then
        case_label="${1}"
        target_directory="${2}"
    else
        echo "ERROR: move_result_to_dir expects a case label and a target directory. Got:"
        echo "  Case label: ${case_label}"
        echo "  target_directory label: ${target_directory}"
        exit 1
    fi

    mv ${ff_solver}.log "${target_directory}/${case_name}_ff.log"
    mv ${pm_solver}.log "${target_directory}/${case_name}_pm.log"

    #ls darcy-iterative
    #find "darcy-iterative*.vtu" -type f -exec ls -al {} \; | sort -nr -k5 | head -n 1
    stokes_final_vtu=$(find .  -maxdepth 1 -iname "stokes-iterative*.vtu" -type f -exec ls {} \; | sort -r | head -n1)
    mv ${stokes_final_vtu} "${target_directory}/stokes-iterative-final.vtu"
    darcy_final_vtu=$(find .  -maxdepth 1 -iname "darcy-iterative*.vtu" -type f -exec ls {} \; | sort -r | head -n1)
    mv ${darcy_final_vtu} "${target_directory}/darcy-iterative-final.vtu"

    for f in $(find .  -maxdepth 1 -iname "precice*-iterations.log" -type f -exec ls {} \+);
    do
        mv ${f} "${target_directory}"
    done

    clean_up_dir
}

##################
# Stokes
##################

run_test_case ./serial-implicit-stokes-first/precice-config.xml
move_result_to_dir "flow-over-box" "serial-implicit-stokes-first/"

run_test_case ./serial-implicit-darcy-first/precice-config.xml
move_result_to_dir "flow-over-box" "serial-implicit-darcy-first/"

run_test_case ./parallel-implicit/precice-config.xml
move_result_to_dir "flow-over-box" "parallel-implicit/"

##################
# Navier-Stokes
##################

run_test_case ./serial-implicit-stokes-first/precice-config.xml
move_result_to_dir "flow-over-box" "serial-implicit-stokes-first/"

run_test_case ./serial-implicit-darcy-first/precice-config.xml
move_result_to_dir "flow-over-box" "serial-implicit-darcy-first/"

run_test_case ./parallel-implicit/precice-config.xml
move_result_to_dir "flow-over-box" "parallel-implicit/"


exit 0