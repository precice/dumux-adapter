#!/usr/bin/env bash

# Runs clang-format on all *.cc and *.hh files
find appl/ -regex '.*\.\(cc\|hh\)' -exec clang-format-10 -style=file -i {} \;