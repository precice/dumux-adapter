# DuMuX-preCICE examples

The source code of the DuMuX-preCICE examples can be found in this directory. After the configure and build step of the adapter, the corresponding build scripts etc. can be found in `build-cmake/examples`. This means for modifying the examples one has to edit the files in the `examples/` directory. For building and running the examples one has to go to the corresponding directory in `build-cmake/examples`.

The examples often have two input files that should be passed to the executables:

1. An `.input` file which is a DuMUX input file describing the simulation setting, e.g., pressure, name of output files or mesh size.
2. One or several `.xml` files which describe preCICE's coupling configuration, e.g. mapping types, data acceleration etc. These files are not provided for monolithic test cases since a preCICE configuration is only needed for partitioned couplings.

**Note:** The examples described here are used for testing the correctness of the adapter as well. These tests are defined in the `CMakeLists.txt` file of the examples and relevant scripts and reference data reside in the `test/` subdirectory.

## Dummy solver

The dummy solver reside in `examples/dummysolver`. The solver does not solve any equations, but uses the DuMuX adapter to communicate some data between two instances of the dummy solver. The exchanged data is also checked within the dummy solvers as running the dummy solver is part of test implemented tests.

The dummy solver can be used as an example on how to use the adapter, but also for developing and debugging a new code.

The executable of the dummy solver is called `dumuxprecice_dummysolver` and one can start two instances via the provided `Allrun.sh` script. Please also refer to the `Allrun.sh` script if you are curious on how to start the solver manually. Several typical preCICE dummy solver parameters, such as the solver name or the mesh name are give as run-time parameters via the DuMuX parameter system. These parameters could also be provided by a DuMuX `.input` file instead of using the command line.
