# DuMuX-preCICE examples

The source code of the DuMuX-preCICE examples can be found in this directory. After the configure and build step of the adapter, the corresponding build scripts etc. can be found in `build-cmake/examples`. This means for modifying the examples one has to edit the files in the `examples/` directory. For building and running the examples one has to go to the corresponding directory in `build-cmake/examples`.

The examples often have two input files that should be passed to the executables:

1. An `.input` file which is a DuMUX input file describing the simulation setting, e.g., pressure, name of output files or mesh size.
2. One or several `.xml` files which describe preCICE's coupling configuration, e.g. mapping types, data acceleration etc. These files are not provided for monolithic test cases since a preCICE configuration is only needed for partitioned couplings.

## Dummy solver

The dummy solver reside in `examples/dummysolver`. The solver does not solve any equations, but uses the DuMuX adapter to communicate some data between two instances of the dummy solver.

The dummy solver is an example on how to use the adapter, but also for developing and debugging a new code.

To start the two dummy participants, run following commands respectively in two terminals from the folder `build-cmake/examples`:

```bash
./dummy_participantOne params_one.input
```

and

```bash
./dummy_participantTwo params_two.input
```
