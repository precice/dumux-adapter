#/usr/bin/env bash

DUNEVERSION="${1}"
DUMUXVERSION="${2}"

if [[ $# -lt 2 ]]; then
    echo "DUNE or DUMUX version missing"
    exit 1
fi

echo "DUNE version ${DUNEVERSION}"
echo "DUMUX version ${DUMUXVERSION}"

for module in 'common' 'geometry' 'grid' 'localfunctions' 'istl' 'subgrid'
do
    echo "Checking out module: ${module}"
    if [[ ${module} == 'subgrid' ]]; then
        git clone --depth 1 https://gitlab.dune-project.org/extensions/dune-${module}.git -b releases/$DUNEVERSION
    else
        git clone --depth 1 https://gitlab.dune-project.org/core/dune-${module}.git -b releases/$DUNEVERSION
    fi
done

if [[ "$DUMUXVERSION" != "none" ]]; then
    echo "Checking out module: dumux"
    if [[ "$DUMUXVERSION" == "master" ]]; then
    git clone --depth 1 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git -b master
    else
    git clone --depth 1 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git -b releases/$DUMUXVERSION
    fi
fi
echo "Printing DUNE modules to screen:"
ls
