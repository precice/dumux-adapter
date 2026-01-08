#!/usr/bin/env sh

set -e -u

rm -rf precice-run/

./dumuxprecice_singlecheckpointtest -preCICE.SolverName SolverOne -preCICE.ConfigFileName test.xml -preCICE.MeshName MeshOne > Solver_One.out 2>&1 &
SOLVER_ONE_ID=$!

./dumuxprecice_singlecheckpointtest -preCICE.SolverName SolverTwo -preCICE.ConfigFileName test.xml -preCICE.MeshName MeshTwo > Solver_Two.out 2>&1 &
SOLVER_TWO_ID=$!

wait ${SOLVER_ONE_ID}
if [ $? -ne 0 ]; then
    exit $?
fi
wait ${SOLVER_TWO_ID}
if [ $? -ne 0 ]; then
    exit $?
fi
