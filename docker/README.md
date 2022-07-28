# Docker image

This directory contains a Docker recipe for building an environment for testing the DuMuX-preCICE adapter. The recipe does **not** install the adapter.

The following main dependencies are installed:

- DUNE
- DuMuX
- preCICE

For more details on what is installed, please check the `Dockerfile`.

The versions of the main dependencies can be manipulated when building the Docker image via the arguments with given default values

- `ARG DUNEVERSION=2.7`: DUNE version
- `ARG DUMUXVERSION=3.4`: DuMuX version
- `ARG PRECICEVERSION=2.2.1`: preCICE version
- `ARG PRECICEUBUNTU=focal`: preCICE target platform (here, Ubuntu Focal Fossa)

If one wants to build an image using DUNE 2.8, DuMuX 3.5 and preCICE 2.4.0 on could call it with the following command:

```text
docker build --build-arg DUNEVERSION=2.8 --build-arg DUMUXVERSION=3.5 --build-arg PRECICEVERSION=2.4.0 -t dune-dumux-precice:2.8-3.5-2.4.0.
```
