#!/usr/bin/env sh

set -e -u

rm -rf precice-run/

./test_partitioned_dummysolver -preCICE.SolverName SolverOne -preCICE.ConfigFileName precice-dummy-solver-config.xml -preCICE.MeshName MeshOne > Solver_One.out 2>&1 &
SOLVER_ONE_ID=$!

./test_partitioned_dummysolver -preCICE.SolverName SolverTwo -preCICE.ConfigFileName precice-dummy-solver-config.xml -preCICE.MeshName MeshTwo > Solver_Two.out 2>&1 &
SOLVER_TWO_ID=$!

wait ${SOLVER_ONE_ID}
if [ $? -ne 0 ]; then
    exit $?
fi
wait ${SOLVER_TWO_ID}
if [ $? -ne 0 ]; then
    exit $?
fi