# DuMuX-preCICE adapter

![build and test develop](https://github.com/precice/dumux-adapter/actions/workflows/build-and-test.yml/badge.svg)
![build and test develop with DuMuX masters](https://github.com/precice/dumux-adapter/actions/workflows/build-and-test-dumux-master.yml/badge.svg)

This repository provides a [DuMuX](https://dumux.org/)-specific adapter to couple to other codes using [preCICE](https://precice.org/). The source code of the adapter was formerly stored [in a repository on the IWS GitLab](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-precice).

## Documentation

Find the [user documentation on the preCICE website](https://precice.org/adapter-dumux.html).

You can also generate Doxygen documentation of the API from the `docs-api/` directory.

## Structure of the repository

Note that this repository is a [DUNE module](https://www.dune-project.org/) and thus some parts of the repository structure are given by the typical DUNE module layout.

- `cmake/`: Contains CMake modules for building the adapter. Under normal circumstances you do not need
- `docs/`: User documentation files.
- `docs-api/`: Doxygen-based API documentation configuration files.
- `docker/`: A Docker recipe that creates a container with DUNE, DuMuX and preCICE. The recipe is mainly used for the automated tests. Check the `README.md` in the subdirectory for more details.
- `dumux-precice/`: The preCICE adapter source code and further code for some of the tests and examples.
- `examples/`: Contains examples on how to couple different domains. Some of the examples are taken from DuMuX or are slightly adapted from DuMuX test cases or tutorials. Please check the `README.md` file in this directory and corresponding subdirectories to find further explanations of the examples.
- `scripts/`: Contains useful scripts to run simulations and for checking the code's formatting.
- `test/`: Contains test cases. The directory also contains several DUNE configuration files (`.opts` files) for configuring the project.
