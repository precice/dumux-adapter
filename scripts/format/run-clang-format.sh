#!/usr/bin/env bash
set -e
# Runs clang-format on all *.cc and *.hh files
DIRECTORIES=("dumux-precice/" "examples/" "test/" "postprocessing/")
for DIRECTORY in "${DIRECTORIES[@]}"; do
    #echo "Check files in ${DIRECTORY}"
    for FILE in $(find "${DIRECTORY}" -regex '.*\.\(cc\|hh\)' ); do
        #echo "$FILE"
        clang-format-10 --dry-run -style=file -Werror -i ${FILE}
    done
    if [ $? -ne 0 ]; then
        exit $?
    fi
done
exit 0