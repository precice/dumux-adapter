#!/usr/bin/env bash

DUNE_VERSION="2.8"
DUMUX_VERSIONS=("3.4" "3.5")
PRECICE_VERSIONS=("develop")

for dumux_version in ${DUMUX_VERSIONS[@]}; do
    for precice_version in ${PRECICE_VERSIONS[@]}; do
        echo "Building dumux-precice:${dumux_version}-${precice_version}"
        docker build --build-arg DUNEVERSION=${DUNE_VERSION} --build-arg DUMUXVERSION=${dumux_version} --build-arg PRECICEVERSION=${precice_version} -t precice/dumux-precice:${dumux_version}-${precice_version} --file dockerfile.slim .
        docker push precice/dumux-precice:${dumux_version}-${precice_version}
        yes | docker image prune
    done
done

PRECICE_VERSION="develop"
echo "Building dune-precice:${dumux_version}-${precice_version}"
docker build --build-arg DUNEVERSION=${DUNE_VERSION} --build-arg PRECICEVERSION=${PRECICE_VERSION} -t precice/dune-precice:${DUNE_VERSION}-${PRECICE_VERSION} --file dockerfile_dune-precice.slim .
docker push precice/dune-precice:${DUNE_VERSION}-${PRECICE_VERSION}
yes | docker image prune
