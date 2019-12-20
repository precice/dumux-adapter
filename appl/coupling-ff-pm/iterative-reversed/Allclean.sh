#! /usr/bin/env bash


ff_solver="test_ff_reversed"
pm_solver="test_pm_reversed"


rm -rf "precice-run/"
rm -f "${ff_solver}.log" "${pm_solver}.log" 
rm -f precice-*.log
rm -f precice-*.json

