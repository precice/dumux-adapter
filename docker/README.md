# Docker image

This directory contains a Docker recipes and helper scripts for building an environment for testing the DuMuX-preCICE adapter. The recipes do **not** install the adapter.

The following main dependencies are installed when using the `dockerfile.large` and `dockerfile.slim` recipe:

- DUNE
- DuMuX
- preCICE

The `dockerfile.base` recipe only installs:

- DUNE
- DuMuX

The "slim" image tries to limit the size of the container by compiling preCICE from source and deactivating many optional features such as language binding, PETSc, Python Actions and MPI.

For more details on what is installed, please check the `Dockerfile`.

The versions of the main dependencies can be manipulated when building the Docker image via the arguments with given default values

- `ARG DUNEVERSION=2.9`: DUNE version. Available in all recipes.
- `ARG DUMUXVERSION=3.7`: DuMuX version. Available in all recipes.
- `ARG PRECICEVERSION=develop`: preCICE version. Available only in "large" and "slim" recipes.
- `ARG PRECICEUBUNTU=focal`: preCICE target platform (here, Ubuntu Focal Fossa). Available only in "large" recipe.

If one wants to build a "slim" image using DUNE 2.9, DuMuX from the current master branch and preCICE 3.0.0 one could call it with the following command:

```text
sudo docker build --build-arg DUNEVERSION=2.9 --build-arg DUMUXVERSION=master --build-arg PRECICEVERSION=3.0.0 -t precice/dumux-precice:master-3.0.0 --file dockerfile.slim .
```

The script `rebuild-docker-images.sh` can be used to rebuild the Docker images locally.
