# DuMuX-preCICE adapter documentation

!!! info
    This documentation homepage is under construction.

This repository provides a [DuMuX](https://dumux.org/)-specific adapter to couple to other codes using [preCICE](https://www.precice.org/). You can find the source code of this adapter [on GitHub](https://github.com/precice/dumux-adapter). The source code of the adapter was formerly stored [in a repository on the IWS GitLab](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-precice). If something is unclear or you would want something to be documented better, please open an [issue](https://github.com/precice/dumux-adapter/issues) or [pull request](https://github.com/precice/dumux-adapter/pulls) and let us know.

## API documentation

The interface of the coupling adapter and also the internal (private) interface are documented using Doxygen. In order to build this documentation you need [Doxygen](https://www.doxygen.nl/index.html) installed. After configuring the project using CMake/`dunecontrol` you can build the documentation via navigating to the `build-cmake` directory and building the `doxygen_dumux-precice` target, i.e.,

```text
cd build-cmake
make doxygen_dumux-precice
```

This generates a HTML documentation which you can view in a browser of your choice. It is stored in `build-cmake/doc/doxygen/index.html`.

## How to cite this code

There is no publication related to this code available yet.

### Publications using dumux-precice

- Jaust A., Weishaupt K., Mehl M., Flemisch B. (2020) Partitioned Coupling Schemes for Free-Flow and Porous-Media Applications with Sharp Interfaces. In: Klöfkorn R., Keilegavlen E., Radu F., Fuhrmann J. (eds) Finite Volumes for Complex Applications IX - Methods, Theoretical Aspects, Examples. FVCA 2020. Springer Proceedings in Mathematics & Statistics, vol 323. Springer, Cham. <https://doi.org/10.1007/978-3-030-43651-3_57>
    - Code can be found at: <https://git.iws.uni-stuttgart.de/dumux-pub/jaust2020a>
