#!/usr/bin/env sh

set -e -u

./dummy_participantTwo params_two.input > Solver_Two.out 2>&1 &
SOLVER_TWO_ID=$!

wait ${SOLVER_TWO_ID}
if [ $? -ne 0 ]; then
    exit $?
fi
