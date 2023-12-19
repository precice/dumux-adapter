# Adapter usage

Please check out the examples in the `examples/` directory to get an idea on how to use the adapter.

To use the adapter in a separate DUNE module, set `dumux-precice` as a (suggested) dependency and link targets against the `dumux-precice` library as done in the `examples/` directory. If `libdumux-precice` is built as a static library, additionally discover and link the `precice::precice` library. To build the adapter library as a dynamic library, use the CMake option `-DBUILD_SHARED_LIBS=ON` to build `dumux-precice` and upstream modules.
