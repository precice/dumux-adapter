#! /usr/bin/env bash

targetRoot="serial-implicit"


pressureDifferences=("1e-9")
permeabilities=("1e-6")
alphaBeaversJoseph=("1.0")
#meshSizes=("20" "40" "80")
meshSizes=("20" "40")
#meshSizes=("40")
hasInertiaTerms=("true" "false")
#hasInertiaTerms=("true")
#preciceRelativeTolerance=("1e-2" "1e-3" "1e-4" "1e-5" "1e-6" "1e-7" "1e-8")
preciceRelativeTolerance=("1e-2" "1e-4" "1e-6" "1e-8")
preciceConfigBase=("base-serial-implicit.xml" "base-serial-implicit-darcy-first.xml")

#solver_input="params.input"
inputTemplate="fvca-iterative-base.input"

ff_solver="test_ff_reversed"
pm_solver="test_pm_reversed"


rm -f *.csv *.pvd *.vtu *.log *.txt *.json

#echo "" > ${configurations}
basedir="${PWD}"
for hasInertiaTerms in "${hasInertiaTerms[@]}"; do
  for preciceBase in "${preciceConfigBase[@]}"; do
    for preciceRelTol in "${preciceRelativeTolerance[@]}"; do
      for dp in "${pressureDifferences[@]}"; do
        for permeability in "${permeabilities[@]}"; do
          for alpha in "${alphaBeaversJoseph[@]}"; do
            for mesh in "${meshSizes[@]}"; do
              i=$((i+1))
              
              # Check if Stokes or Navier-Stokes        
              flowProblemName="stokes"
              echo ${hasInertiaTerms}
              if [[ "${hasInertiaTerms}" == "true" ]]; then
                flowProblemName="navier-stokes"
              fi
            
              # Generate name of test case and create directories
              casename="${flowProblemName}-${preciceRelTol}-${mesh}-${alpha}-${permeability}-${dp}"
              if [[ ${preciceBase} == *"darcy-first"* ]]; then
                casename="${casename}-darcy-first"
              fi
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
                  
              preciceXML="${casename}.xml"
              echo "${preciceBase}"
              sed -e "s/RELTOL/${preciceRelTol}/g" \
                  "${preciceBase}" > ${preciceXML}                  
                  
              rm -rf "precice-run/"
#              ff_cmd="./${ff_solver} - ${preciceXML}"
              echo "${ff_cmd}"
              time ./${ff_solver} ${inputFile} - ${preciceXML} > ${ff_solver}.log 2>&1 &
              PIDFluid=$!
              time ./${pm_solver} ${inputFile} - ${preciceXML} > ${pm_solver}.log 2>&1 &
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


              targetDir="${targetRoot}/${casename}"
              mkdir -p "${targetDir}"
              mv *.csv ${targetDir}/
              mv *.vtu ${targetDir}/
              mv *.pvd ${targetDir}/
              mv *.log ${targetDir}/
              mv *.txt ${targetDir}/
              mv *.json ${targetDir}/
              mv ${preciceXML} ${targetDir}/
              mv "${inputFile}" ${targetDir}/ 
              rm -f *.csv *.pvd *.vtu *.log *.txt *.json
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
done

echo "In total ${i} cases were run"





