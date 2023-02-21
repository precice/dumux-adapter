# Structure of the repository

Note that the DuMuX-preCICE adapter is a [DUNE module](https://www.dune-project.org/) and thus some parts of the repository structure are given by the typical DUNE module layout.

- `cmake/`: Contains CMake modules for building the adapter. Under normal circumstances you do not need
- `examples/`: Contains examples on how to couple different domains. Some of the examples are taken from DuMuX or are slightly adapted from DuMuX test cases or tutorials. Please check the `README.md` file in this directory and corresponding subdirectories to find further explanations of the examples. Additional examples can be found in the `test/` directory.
- `doc/`: Additional documentation.
- `docker/`: A Docker recipe that creates a container with DUNE, DuMuX and preCICE. The recipe is mainly used for the automated tests. Check the `README.md` in the subdirectory for more details.
- `dumux-precice/`: The preCICE adapter source code and further code for some of the tests and examples.
- `scripts/`: Contains useful scripts to run simulations and for checking the code's formatting.
- `test/`: Contains test cases and reference solutions (`reference-solutions/`). The directory also contains several DUNE configuration files (`.opts` files) for configuring the project.
