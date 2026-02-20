---
title: Use the DuMuX adapter
permalink: adapter-dumux-use.html
keywords: DuMuX, DUNE, C++
summary: How to use the DuMuX adapter for building your own coupled solver.
---

To understand how the adapter is used, see the tutorial [free-flow-over-porous-media](https://github.com/precice/tutorials/tree/develop/free-flow-over-porous-media). Additionally, the solver [macro-dumux](https://github.com/precice/tutorials/tree/develop/two-scale-heat-conduction/macro-dumux) uses the adapter to couple a heat conduction problem to micro-scale simulations.

### Configuration

The adapter is configured via the group called `[precice-adapter-config]` in the DuMux runtime parameter file with the default name `params.input`. The input files for the dummy simulation under `examples/dummysolver` provided examples of the parameters to be configured. This configuration follows the nomenclature of the [preCICE adapter configuration schema](https://github.com/precice/preeco-orga/tree/main/adapter-config-schema). It does not adhere to the schema completely because the configuration is done via a `.input` file instead of a JSON or YAML file.

### Use the API

#### Set coupling mesh

To inform preCICE the coupling mesh, vertices coordinates and `meshName` are to be set with `setSurfaceMesh()` or `setVolumeMesh()`.

Both functions require:

- `meshName`, which is the name of the mesh to be set
- `coupledDumuxIDs`, which can be indices of the finite-volume face in a surface coupling case or indices of elements in a volume coupling case, and
- `positions`,which are the spatial coordinates that the previous IDs represent.

In addition, `setVolumeMesh()` requires `gridView`, which is used to filter the `overlap` and `ghost` cells in a distributed setting. This is meaningful for defining the correct number of elements that are coupled.

In a distributed surface coupling case, there can also be `overlap` and `ghost` cells on the coupling interface. The filtering of these indices are not supported in the adapter yet, as this also requires information of the geometry. It's required, for user, to only add cell faces belonging to `interior` elements to the `coupledDumuxIDs` before calling `setSurfaceMesh()`.

#### Exchange data

In a parallel setting, the values on the `overlap` and `ghost` elements need to be updated manually via subdomain communications, e.g. `communicate()` provided by DUNE-Grid, after the values of the coupling variables have been set by data from another solver. This is needed since only `interior` elements or cell faces attached to `interior` cells are included in coupled mesh.

### Additional build guideline

To use the adapter in a separate DUNE module, call `dune_enable_all_packages()` in the root `CMakeLists.txt` of the application module. If `libdumux-precice` is built as a static library, preCICE needs to be explicitly discovered with `find_package` as done in the root `CMakeLists.txt` of the adapter. To build the adapter library as a dynamic library, use the CMake option `-DBUILD_SHARED_LIBS=ON` to build `dumux-precice` and upstream modules.
