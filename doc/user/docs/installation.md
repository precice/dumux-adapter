# Installation

The DuMuX-preCICE adapter is a DUNE module named `dumux-precice` which can be build using the [DUNE build system](https://www.dune-project.org/doc/installation/). The DUNE build system is build on top of CMake and comes with various tools that make installation and management of DUNE modules easier. Therefore, it is recommended to install `dumux-precice` using `dunecontrol`. Please check out the [DUNE installation instructions](https://www.dune-project.org/doc/installation/) to get an overview over the `dunecontrol` tools and how DUNE modules work.

## Prerequisites

- DuMuX **newer** than 3.2

    - Builds using the current `master` branch of DuMuX might fail.
    - If you run into trouble with a new DuMuX release, please open an issue in the repository and add the error message that you receive.

- preCICE >=3.0.0

    - The adapter is build via the DUNE build system that is based on CMake. Thus, the CMake [link instructions for preCICE](https://precice.org/installation-linking.html#cmake) apply.

- `wget` or `git` to download the DuMuX-preCICE adapter.
- Optional: [`dune-subgrid`](https://www.dune-project.org/modules/dune-subgrid/) allows for modified grid geometries.

The DuMuX-preCICE adapter should build fine if DuMuX, preCICE and their dependencies are installed.

## Detailed installation steps

1. Install [DuMuX](https://dumux.org/) and the needed dependencies. The easiest way is to follow [DuMuX's installation instructions](https://dumux.org/docs/doxygen/master/installation.html). The DuMuX project provides a script that installs and DuMuX and the DUNE modules required by DuMuX. This means, after installing DuMuX via the provided script you should be good to go to use the DuMuX-preCICE adapter.

    After the installation you should have a root directory that contains the base DUNE modules, i.e. a  number of directories named like `dune-common`, `dune-geometry` etc., and a directory called `dumux`.

    - Note that extended features of DuMuX or the DuMuX-preCICE adapter may need additional DUNE modules.

2. Download the DuMuX-preCICE adapter to the same directory as the DUNE modules and the `dumux` folder. At the moment there are no adapter releases, besides an outdated `v0.1` release, such the best way is to checkout the `develop` branch of the adapter.

    ```text
    git clone -b develop https://github.com/precice/dumux-adapter.git
    ```

    You can also try to clone the repository via SSH:

    ```text
    git clone -b develop git@github.com:precice/dumux-adapter.git
    ```

3. Verify that the `dumux-adapter` folder is in the same directory as the DUNE module folders and the `dumux` folder.

4. Build and configure the adapter using `dunecontrol`. While being in the directory mentioned in the previous step via calling

    ```text
    dune-common/bin/dunecontrol --only=dumux-precice all
    ```

    After the build and configure step a new directory `build-cmake` was created inside the `dumux-precice` directory.

    You can configure the build and configuration process using advanced options by manipulating CMake variables. `dunecontrol` allows to pass an options file for that

    ```bash
    dune-common/bin/dunecontrol --opts=OPTSFILE.opts --only=dumux-precice all
    ```

    There is an `opts`-file provided by the adapter that resides in `test/`. You can use it as

    ```bash
    dune-common/bin/dunecontrol --opts=dumux-adapter/test/cmake-test.opts --only=dumux-precice all
    ```

    This provided `cmake-test.opts` file turns off some system-dependent optimizations such that the tests create comparable results on different computers.

    To use the adapter in a separate DUNE module, we recommend building the adapter as a shared library. To do so, use the CMake option `-DBUILD_SHARED_LIBS=ON` to build the adapter and upstream modules. The DuMux install script uses `dumux/cmake.opts`, which already sets this option.
    Note that to change this setting it may be required to clear the CMake caches in `build-cmake/CMakeCache.txt`.

    For more ways do manipulate/adapt the build and configuration step, please consult the `dunecontrol` documentation.

5. Optional, but recommended: Build all tests to verify the installation. For this navigate in the `build-cmake/` directory and build the `build_tests` target.

    ```bash
    cd dumux-adapter/build-cmake
    make -j1 build_tests
    ```

    You may speed up the build process by using more than one build job, e.g., use `make -j4` in order to build with for processes at the same time.

    Afterwards you can run the tests from the `build-cmake` directory using

    ```bash
    ctest
    ```

    If any tests fails, you should verify if something went wrong with the installation.

There are advanced ways of managing DUNE modules, e.g. using the environment variable `DUNE_CONTROL_PATH`, that are beyond the scope of this short documentation. You can find more information in the [DUNE FAQ](https://www.dune-project.org/doc/installation/#faq).
