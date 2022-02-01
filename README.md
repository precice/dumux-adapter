# DuMuX-preCICE coupling

This repository provides a [DuMuX](https://dumux.org/)-specific adapter to couple to other codes using [preCICE](https://www.precice.org/). You can find the source code of this adapter [on the IWS GitLab](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-precice).

## Prerequisites

- DuMuX **newer** than 3.2
    - Builds using the current `master` branch of DuMuX might fail.
- preCICE >=2.0.0
    - The adapter is build via the DUNE build system that is based on CMake. Thus, the CMake [link instructions for preCICE](https://precice.org/installation-linking.html#cmake) apply.

## Structure of the repository

- `cmake/`: Contains CMake modules for building the adapter. Under normal circumstances you do not need
- `examples/`: Contains examples on how to couple different domains. Some of the examples are taken from DuMuX or are slightly adapted from DuMuX test cases or tutorials. Please check the `README.md` file in this directory and corresponding subdirectories to find further explanations of the examples. Additional examples can be found in the `test/` directory.
- `doc/`: Additional documentation.
- `dumux-precice/`: The preCICE adapter source code and further code for some of the tests and examples.
- `scripts/`: Contains useful scripts to run simulations and for checking the code's formatting.
- `test/`: Contains test cases and reference solutions (`reference-solutions/`). Test cases as divided into `monolithic/` and `partitioned` (=preCICE) test cases. The directory also contains several DUNE configuration files (`.opts` files) for configuring the project.

## Running tests

1. Configure the DUNE module

    ```bash
    dunecontrol --opts=OPTSFILE.opts --only=dumux-precice all
    ```

    You may use one of the `opts`-file provided by the adapter that resides in `test/`, e.g.,

    ```bash
    dunecontrol --opts=dumux-precice/test/cmake-test.opts --only=dumux-precice all
    ```

2. Build all test. For this navigate in the `build-cmake/` directory and build the `build_tests` target.

    ```bash
    cd dumux-precice/build-cmake
    make -j1 build_tests
    ```

    You may speed up the build process by using more than one build job, e.g., use `make -j4` in order to build with for processes at the same time.

3. Run the tests from the `build-cmake` directory

    ```bash
    ctest
    ```

## Documentation

### User documentation

At the moment the documentation is provided by the API documentation (see below) as well as the test and example cases. If something is unclear or you would want something to be documented better, please open an [issue](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-precice/-/issues) and let us know.

### API documentation

The interface of the coupling adapter and also the internal (private) interface are documented using Doxygen. In order to build this documentation you need [Doxygen](https://www.doxygen.nl/index.html) installed. After configuring the project using CMake/`dunecontrol` you can build the documentation via navigating to the `build-cmake` directory and building the `doxygen_dumux-precice` target, i.e.,

```text
cd build-cmake
make doxygen_dumux-precice
```

This generates a HTML documentation which you can view in a browser of your choice. It is stored in `build-cmake/doc/doxygen/index.html`.

## Publications

### How to cite this code?

There is no publication related to this code available yet.

### Publications using dumux-precice

- Jaust A., Weishaupt K., Mehl M., Flemisch B. (2020) Partitioned Coupling Schemes for Free-Flow and Porous-Media Applications with Sharp Interfaces. In: Klöfkorn R., Keilegavlen E., Radu F., Fuhrmann J. (eds) Finite Volumes for Complex Applications IX - Methods, Theoretical Aspects, Examples. FVCA 2020. Springer Proceedings in Mathematics & Statistics, vol 323. Springer, Cham. <https://doi.org/10.1007/978-3-030-43651-3_57>
    - Code can be found at: https://git.iws.uni-stuttgart.de/dumux-pub/jaust2020a
