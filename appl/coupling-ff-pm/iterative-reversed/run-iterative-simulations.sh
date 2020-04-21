#! /usr/bin/env bash


pressureDifferences=("1e-9")
permeabilities=("1e-6")
alphaBeaversJoseph=("1.0")
meshSizes=("20" "40" "80")
hasInertiaTerms=("true" "false")
precice_config=("precice-config-serial-implicit-reversed.xml" "precice-config-serial-implicit-reversed-darcy-first.xml")

#solver_input="params.input"
inputTemplate="fvca-iterative-base.input"

ff_solver="test_ff_reversed"
pm_solver="test_pm_reversed"


rm -f *.csv *.pvd *.vtu *.log *.txt *.json

#echo "" > ${configurations}
basedir="${PWD}"
for preciceFile in "${precice_config[@]}"; do
  for hasInertiaTerms in "${hasInertiaTerms[@]}"; do
    for dp in "${pressureDifferences[@]}"; do
      for permeability in "${permeabilities[@]}"; do
        for alpha in "${alphaBeaversJoseph[@]}"; do
          for mesh in "${meshSizes[@]}"; do
            i=$((i+1))
            
            # Check if Stokes or Navier-Stokes        
            flowProblemName="stokes"
  #          echo ${hasInertiaTerms}
            if [[ "${hasInertiaTerms}" == "true" ]]; then
              flowProblemName="navier-stokes"
            fi
          
            # Generate name of test case and create directories
            casename="${flowProblemName}-${mesh}-${alpha}-${permeability}-${dp}"
            echo "${casename}"

            # Setting up input file          
            inputFile="${casename}.input"
            sed -e s/MESHSIZE/"${mesh}"/g \
                -e "s/FLOWPROBLEMNAME/${flowProblemName}/g" \
                -e "s/PRESSUREDIFF/${dp}/g" \
                -e "s/PERM/${permeability}/g" \
                -e "s/ALPHA/${alpha}/g" \
                -e "s/HASINERTIATERMS/${hasInertiaTerms}/g" \
                -e "s/CASENAME/${casename}/g" \
                "${inputTemplate}" > ${inputFile}
                
            # Running simulation
            #solverCmd="./${solver} ${inputFile}"
            #echo "$solverCmd"
   
  #           $(${solverCmd} > solver.log)                   
#            $(time ./${solver} "${inputFile}" > solver.log)
  #          $(./${solverCmd} ${inputFile} > solver.log)


            rm -rf "precice-run/"
            ff_cmd="./${ff_solver} - ${precice_config}"
            echo "${ff_cmd}"
            ./${ff_solver} - ${preciceFile} > ${ff_solver}.log 2>&1 &
            PIDFluid=$!
            ./${pm_solver} - ${preciceFile} > ${pm_solver}.log 2>&1 &
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


            mkdir -p "${casename}"
            mv *.csv ${casename}
            mv *.vtu ${casename}
            mv *.pvd ${casename}
            mv *.log ${casename}
            mv *.txt ${casename}
            mv *.json ${casename}
            cp ${preciceFile} ${casename}
            mv "${inputFile}" ${casename} 
            #cd ${casename}
  #          ln -s "../${solver}" "${solver}"
            #cp "../${solver}" .

            # Go back to root dir
            cd "${basedir}"
          done
        done
      done
    done
  done
done

echo "In total ${i} cases were run"





