# Example cases

The source code of the DuMuX-preCICE examples can be found in the subdirectory `examples/`. After the configure and build step of the adapter, the corresponding build scripts etc. can be found in `build-cmake/examples`. This means for modifying the examples one has to edit the files in the `examples/` directory. For building and running the examples one has to go to the corresponding directory in `build-cmake/examples`.

The examples often have two input files that should be passed to the executables:

1. An `.input` file which is a DuMUX input file describing the simulation setting, e.g., pressure, name of output files or mesh size.
2. One or several `.xml` files which describe preCICE's coupling configuration, e.g. mapping types, data acceleration etc. These files are not provided for monolithic test cases since a preCICE configuration is only needed for partitioned couplings.

!!! info "Test cases"
    The examples described here are used for testing the correctness of the adapter as well. These tests are defined in the `CMakeLists.txt` file of the examples and relevant scripts and reference data reside in the `test/` subdirectory.

## Free and porous medium flow

All examples of coupling free and porous-medium flow reside in `examples/ff-pm` and corresponding subdirectories in there. We use the abbreviation `ff` for things concerning the free flow and `pm` for things concerning the porous-medium part.

### Flow over a porous-medium

We implemented a simple test case with a pressure driven flow in a free-flow subdomain (top of domain). At the bottom of the free-flow subdomain the subdomain is connected to the porous-medium subdomain. All other boundaries of the porous-medium subdomain are walls (no-flow boundaries).

|  | 2D | 3D |
| --- | --- | --- |
| Directory | `examples/flow-over-square-2d` | `examples/flow-over-cube-3d` |
| Name of executable(s) | `ff_flow_over_square_2d` and `pm_flow_over_square_2d`| `ff_flow_over_cube_3d` and `pm_flow_over_cube_3d` |

The 2D test case comes with three preCICE configurations for parallel-implicit coupling (`pi`), serial-implicit coupling (`si`) with either running the flow simulation first (`free-flow-first`) or second (`free-flow-second`). The 3D test case comes with only one preCICE configuration file.

## Dummy solver

The dummy solver reside in `examples/dummysolver`. The solver does not solve any equations, but uses the DuMuX adapter to communicate some data between two instances of the dummy solver. The exchanged data is also checked within the dummy solvers as running the dummy solver is part of test implemented tests.

The dummy solver can be used as an example on how to use the adapter, but also for developing and debugging a new code.

The executable of the dummy solver is called `dumuxprecice_dummysolver` and one can start two instances via the provided `Allrun.sh` script. Please also refer to the `Allrun.sh` script if you are curious on how to start the solver manually. Several typical preCICE dummy solver parameters, such as the solver name or the mesh name are give as run-time parameters via the DuMuX parameter system. These parameters could also be provided by a DuMuX `.input` file instead of using the command line.
