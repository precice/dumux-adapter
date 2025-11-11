#!/usr/bin/env sh

set -e -u

rm -rf precice-run/

./adapterConfigTest paramsSolverOne.input > Solver_One.out 2>&1 &
SOLVER_ONE_ID=$!

./adapterConfigTest paramsSolverTwo.input > Solver_Two.out 2>&1 &
SOLVER_TWO_ID=$!

wait ${SOLVER_ONE_ID}
if [ $? -ne 0 ]; then
    exit $?
fi
wait ${SOLVER_TWO_ID}
if [ $? -ne 0 ]; then
    exit $?
fi
