#!/usr/bin/env bash

DUNE_VERSION="2.8"
DUMUX_VERSIONS=("3.4" "3.5")
PRECICE_VERSIONS=("2.4.0" "2.3.0")

for dumux_version in ${DUMUX_VERSIONS[@]}; do
    for precice_version in ${PRECICE_VERSIONS[@]}; do
        echo "Building dumux-precice:${dumux_version}-${precice_version}"
        docker build --build-arg DUNEVERSION=${DUNE_VERSION} --build-arg DUMUXVERSION=${dumux_version} --build-arg PRECICEVERSION=${precice_version} -t ajaust/dumux-precice:${dumux_version}-${precice_version} --file dockerfile.slim .
        docker push ajaust/dumux-precice:${dumux_version}-${precice_version}
        y | docker image prune
    done
done

PRECICE_VERSION="2.4.0"
echo "Building dune-precice:${dumux_version}-${precice_version}"
docker build --build-arg DUNEVERSION=${DUNE_VERSION} --build-arg PRECICEVERSION=${PRECICE_VERSION} -t ajaust/dune-precice:2.8-${PRECICE_VERSION} --file dockerfile_dune-precice.slim .
docker push ajaust/dune-precice:${dumux_version}-${precice_version}
y | docker image prune