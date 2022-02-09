# DuMuX-preCICE change log

## Not released yet

- 2022-02-09: Made sure all private member of the adapter are suffixed with an underscore.
- 2022-02-01: Add some extra information on the documentation in the `README.md`. Removed old/out-of-date mkdocs documentation from `doc/mkdocs`.
- 2022-01-31: Increased robustness of test scripts.
- 2022-01-31: We now use `diff -w` to compare preCICE's output files for regression tests. In preCICE 2.3.0 the white spaces used have changed which broke our regressions tests.
- 2022-01-26: Renamed `dumupreciceindexwrapper.[hh|cc]` to `dumupreciceindexmapper.[hh|cc]` to be consistent with the class name.
- 2022-01-26: Add and configure Doxygen code documentation of coupling adapter.
- 2022-01-25: Fix code formatting configuration to be close to the original DuMuX code formatting configuration.
- 2022-01-25: Added [description templates](https://docs.gitlab.com/ee/user/project/description_templates.html) for merge requests and issues.
- 2022-01-12: The repository has been restructured. The main changes are:

  - The adapter is now called `CouplingAdapter` and resides in `dumux-precice/`. The build process has been adapted accordingly.
  - Tests case reside in `test/` directory and there in the corresponding subdirectory depending on whether it is a `monolithic`ly or `partitioned`ly coupled test case.
  - Other example cases reside in the directory `examples/`. This is mainly the directory called `appl/` before, but with a new folder structure.
  - The configuration of tests has been changed such that it is possible to build all tests using the `build_tests` target.

  For details check out the merge request [!18 Restructure repository and tests](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-precice/-/merge_requests/18)
- 2022-01-10: Add license file. The code is licensed under GPLv3 without template exception.
- 2022-01-10: Tests run by the CI on DuMuX `master` are allowed to fail.
- 2022-01-10: Added `CHANGELOG.md` to track changes of the adapter.

## v0.1

This marks the initial release of the DuMuX-preCICE adapter.

- Should represent the state of adapter used in publication [Jaust2020a](https://git.iws.uni-stuttgart.de/dumux-pub/jaust2020a).
- Requires preCICE 1.6