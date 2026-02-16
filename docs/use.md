---
title: Use the DuMuX adapter
permalink: adapter-dumux-use.html
keywords: DuMuX, DUNE, C++
summary: How to use the DuMuX adapter for building your own coupled solver.
---

To understand how the adapter is used, see the tutorial [free-flow-over-porous-media](https://github.com/precice/tutorials/tree/develop/free-flow-over-porous-media). Additionally, the solver [macro-dumux](https://github.com/precice/tutorials/tree/develop/two-scale-heat-conduction/macro-dumux) uses the adapter to couple a heat conduction problem to micro-scale simulations.

The adapter is configured via the group called `[precice-adapter-config]` in the DuMux runtime parameter file with the default name `params.input`. The input files for the dummy simulation under `examples/dummysolver` provided examples of the parameters to be configured. This configuration follows the [preCICE adapter configuration schema](https://github.com/precice/preeco-orga/tree/main/adapter-config-schema).

To use the adapter in a separate DUNE module, call `dune_enable_all_packages()` in the root `CMakeLists.txt` of the application module. If `libdumux-precice` is built as a static library, preCICE needs to be explicitly discovered with `find_package` as done in the root `CMakeLists.txt` of the adapter. To build the adapter library as a dynamic library, use the CMake option `-DBUILD_SHARED_LIBS=ON` to build `dumux-precice` and upstream modules.
