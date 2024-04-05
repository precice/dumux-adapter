# API documentation

The interface of the coupling adapter and also the internal (private) interface are documented using Doxygen. In order to build this documentation you need [Doxygen](https://www.doxygen.nl/index.html) installed. After configuring the project using CMake/`dunecontrol` you can build the documentation via navigating to the `build-cmake` directory and building the `doxygen_dumux-precice` target, i.e.,

```text
cd build-cmake
make doxygen_dumux-precice
```

This generates a HTML documentation which you can view in a browser of your choice. It is stored in `build-cmake/doc/doxygen/index.html`.
