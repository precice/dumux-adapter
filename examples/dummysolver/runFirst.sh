#!/usr/bin/env sh

set -e -u

rm -rf precice-run/

./dummy_participantOne params_one.input > Solver_One.out 2>&1 &
SOLVER_ONE_ID=$!

wait ${SOLVER_ONE_ID}
if [ $? -ne 0 ]; then
    exit $?
fi
