#!/usr/bin/env bash

path_to_clean=$1
echo ${path_to_clean}

rm -f ${path_to_clean}/*.vtu
rm -f ${path_to_clean}/*.pvd

exit 0