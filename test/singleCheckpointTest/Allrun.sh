#!/usr/bin/env sh

set -e -u

rm -rf precice-run/

./dumuxprecice_singlecheckpointtest params_one.input > Solver_One.out 2>&1 &
SOLVER_ONE_ID=$!

./dumuxprecice_singlecheckpointtest params_two.input > Solver_Two.out 2>&1 &
SOLVER_TWO_ID=$!

wait ${SOLVER_ONE_ID}
if [ $? -ne 0 ]; then
    exit $?
fi
wait ${SOLVER_TWO_ID}
if [ $? -ne 0 ]; then
    exit $?
fi
