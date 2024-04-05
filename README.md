# DuMuX-preCICE adapter

![build and test develop](https://github.com/precice/dumux-adapter/actions/workflows/build-and-test.yml/badge.svg)
![build and test develop with DuMuX masters](https://github.com/precice/dumux-adapter/actions/workflows/build-and-test-dumux-master.yml/badge.svg)

This repository provides a [DuMuX](https://dumux.org/)-specific adapter to couple to other codes using [preCICE](https://www.precice.org/). The source code of the adapter was formerly stored [in a repository on the IWS GitLab](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-precice).

## Structure of the repository

Note that this repository is a [DUNE module](https://www.dune-project.org/) and thus some parts of the repository structure are given by the typical DUNE module layout.

- `cmake/`: Contains CMake modules for building the adapter. Under normal circumstances you do not need
- `examples/`: Contains examples on how to couple different domains. Some of the examples are taken from DuMuX or are slightly adapted from DuMuX test cases or tutorials. Please check the `README.md` file in this directory and corresponding subdirectories to find further explanations of the examples. Additional examples can be found in the `test/` directory.
- `doc/`: Additional documentation.
- `docker/`: A Docker recipe that creates a container with DUNE, DuMuX and preCICE. The recipe is mainly used for the automated tests. Check the `README.md` in the subdirectory for more details.
- `dumux-precice/`: The preCICE adapter source code and further code for some of the tests and examples.
- `scripts/`: Contains useful scripts to run simulations and for checking the code's formatting.
- `test/`: Contains test cases and reference solutions (`reference-solutions/`). The directory also contains several DUNE configuration files (`.opts` files) for configuring the project.

## Documentation

### User documentation

The main user documentation is currently not available on an online service, but can be found in the `docs/user/docs` directory. Additionally, one may find interesting information in the API documentation (see below) as well as the test and example cases that are provided with this repository. If something is unclear or you would want something to be documented better, please open an [issue](https://github.com/precice/dumux-adapter/issues) and let us know.

### API documentation

The interface of the coupling adapter and also the internal (private) interface are documented using Doxygen. In order to build this documentation you need [Doxygen](https://www.doxygen.nl/index.html) installed. After configuring the project using CMake/`dunecontrol` you can build the documentation via navigating to the `build-cmake` directory and building the `doxygen_dumux-precice` target, i.e.,

```text
cd build-cmake
make doxygen_dumux-precice
```

This generates a HTML documentation which you can view in a browser of your choice. It is stored in `build-cmake/doc/doxygen/index.html`.
