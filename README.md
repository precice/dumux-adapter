# DuMuX-preCICE coupling

This repository provides a DuMuX specific adapter to couple to other codes using [preCICE](https://www.precice.org/)
## Prerequisites

- DuMuX **newer** than 3.2
    - Some of
- preCICE 1.6.X
    - We are working on supporting newer versions of preCICE at the moment.

## Running partitioned coupling

At the moment the following has to be fulfilled

1. `test_freeflow` and `test_solidenergy` have to be in the same directory
1. A preCICE configuration file named `precice-config.xml` must be present in the same directory as the two solvers.
1. Example configurations can be found in `/build-cmake/appl/conjugateheattransfer/iterative/precice-configs`. The config files have to be copied and renamed.

**Note**: There is no guarantee that the config file work at the moment. (2019-04-10)

## Notes regarding monolithic examples

We provide some monolithic simulations for comparison. Most of these examples are based on a DuMuX branch that "reverse" the information exchange direction of the monolithic coupling. This means, it allows to prescribe pressure on the coupling boundary of the porous medium and the velocity/mass flux on the free flow. The *standard* approach of DuMuX is the opposite coupling direction, i.e. velocity/mass flux prescriped on the coupling boundary of the porous medium and pressure prescribed on the coupling boundary of the free flow. It is not clear whether the *non-standard* approach can/will be merged and maintained in the DuMuX main branch.

As a consequence, the provided monolithic examples may not work. The partitioned/iterative examples should work after all.
## Publications

This list contains publications that are using this repository:

- Jaust A., Weishaupt K., Mehl M., Flemisch B. (2020) Partitioned Coupling Schemes for Free-Flow and Porous-Media Applications with Sharp Interfaces. In: Klöfkorn R., Keilegavlen E., Radu F., Fuhrmann J. (eds) Finite Volumes for Complex Applications IX - Methods, Theoretical Aspects, Examples. FVCA 2020. Springer Proceedings in Mathematics & Statistics, vol 323. Springer, Cham. https://doi.org/10.1007/978-3-030-43651-3_57
    - Code can be found at: https://git.iws.uni-stuttgart.de/dumux-pub/jaust2020a

## Further reading

There is more documentation in `doc/mkdocs`.
