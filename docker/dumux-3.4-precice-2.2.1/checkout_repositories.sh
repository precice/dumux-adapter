#/usr/bin/env sh

DUNEVERSION="${1}"
DUMUXVERSION="${2}"

echo "DUNE version ${DUNEVERSION}"
echo "DUMUX version ${DUMUXVERSION}"

for module in ('common', 'geometry', 'grid', 'localfunctions', 'istl', 'subgrid')
do
    echo "Checking out module: ${module}"
    if [[ module == 'subgrid']]; then
        git clone https://git.imp.fu-berlin.de/agnumpde/dune-${module}.git -b releases/$DUNEVERSION
    else
        git clone https://gitlab.dune-project.org/core/dune-${module}.git -b releases/$DUNEVERSION
    fi
done

git clone https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git -b releases/$DUMUXVERSION
