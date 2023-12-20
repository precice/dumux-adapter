# Adapter usage

Please check out the examples in the `examples/` directory to get an idea on how to use the adapter.

To use the adapter in a separate DUNE module, call `dune_enable_all_packages()` in `CMakeLists.txt`. If `libdumux-precice` is built as a static library, additionally discover preCICE with `find_package` as done in the root `CMakeLists.txt`. To build the adapter library as a dynamic library, use the CMake option `-DBUILD_SHARED_LIBS=ON` to build `dumux-precice` and upstream modules.
