# DuMuX-preCICE coupling 


## Prerequisites

- DuMuX **newer** than 3.0
- preCICE 1.3
    - Newer versions of preCICE are likely to work, but are not tested at the moment. 

## Running partitioned coupling

At the moment the following has to be fulfilled

1. `test_freeflow` and `test_solidenergy` have to be in the same directory
1. A preCICE configuration file named `precice-config.xml` must be present in the same directory as the two solvers.
1. Example configurations can be found in `/build-cmake/appl/conjugateheattransfer/iterative/precice-configs`. The config files have to be copied and renamed.

**Note**: There is no guarantee that the config file work at the moment. (2019-04-10)

## Further reading

There is more documentation in `doc/mkdocs`.
