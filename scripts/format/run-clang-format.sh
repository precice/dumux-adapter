#!/usr/bin/env sh
set -e
# Runs clang-format on all *.cc and *.hh files
find dumux-precice/ -regex '.*\.\(cc\|hh\)' -exec clang-format-10 --dry-run -style=file -i {} \;
find examples/ -regex '.*\.\(cc\|hh\)' -exec clang-format-10 --dry-run -style=file -i {} \;
find postprocessing/ -regex '.*\.\(cc\|hh\|hpp\|cpp\)' -exec clang-format-10 --dry-run -style=file -i {} \;
