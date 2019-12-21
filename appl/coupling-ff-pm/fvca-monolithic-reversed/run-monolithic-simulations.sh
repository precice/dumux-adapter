#! /usr/bin/env bash

#testcase="VerticalFlow"

pressureDifferences=("1e-9")
permeabilities=("1e-6")
alphaBeaversJoseph=("1.0")
meshSizes=("20" "40" "80")
hasInertiaTerms=("true" "false")

#interpolation=("nearest-neighbor" "nearest-projection")




i=0

#mkdir -p ${testcase}
#cd ${testcase}

solver="fvca-monolithic-reversed"
inputTemplate="fvca-monolithic-base.input"

#cp ../${geomtemplate} .

#configurations="configs-${testcase}.txt"

# Initial clean up

rm -f *.csv *.pvd *.vtu *.log *.txt

#echo "" > ${configurations}
basedir="${PWD}"
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
          solverCmd="./${solver} ${inputFile}"
          #echo "$solverCmd"
 
#           $(${solverCmd} > solver.log)                   
          $(time ./${solver} "${inputFile}" > solver.log)
#          $(./${solverCmd} ${inputFile} > solver.log)


          mkdir -p "${casename}"
          mv *.csv ${casename}
          mv *.vtu ${casename}
          mv *.pvd ${casename}
          mv *.log ${casename}
          mv *.txt ${casename}
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

echo "In total ${i} cases were run"
