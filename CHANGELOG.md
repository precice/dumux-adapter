# DuMuX-preCICE change log

## v2.0.0

- 2023-08-10: Updated the adapter to align with preCICE V3. Updated the examples accordingly.

## v1.0.0

- 2022-09-14: The solver dummy has been cleaned up.
- 2022-09-14: Updated CI to use containers from `precice` namespace on DockerHub. Also makes sure that `dune-precice` containers are build for tests against DuMuX's master branch.
- 2022-09-14: Updated and fixed GitHub to build Docker containers.
- 2022-08-12: Update CI config to rebuild Docker containers if Docker recips have been updated or when the CI config for building the containers has been updated.
- 2022-08-12: Update examples to work with DuMuX releases newer than DuMuX 3.5.
- 2022-08-12: Update CI to test with preCICE 2.5.0.
- 2022-08-11: Updated documentation and removed it from readthedocs.
- 2022-08-09: Remove GitLab-specific files like CI configuration.
- 2022-08-09: Update links to point to the new GitHub repository of the code.
- 2022-08-09: Seperate workflow such that canary builds (builds using DuMuX master) are done independently of the normal CI tests.
- 2022-08-05: Moved CI to GitHub and updated CI workflows. Docker images now use the root user default user. Several Docker images have been added.
- 2022-08-04: Restructure repository to have a more logical directory structure. Code examples live in `examples/`. These codes are also used as tests. Reference data resides in `test/`.
- 2022-08-04: Update CI to also run tests against DuMuX master. Theses tests may fail.
- 2022-08-04: Remove leftover references to monolithic coupling from code.
- 2022-08-02: Fix CMake macro that creates files links to preCICE configuration files.
- 2022-07-29: Updated tests to work with DuMuX 3.5. CI was updated to newer versions of DuMuX and preCICE accordingly.
- 2022-07-28: Remove monolithic test cases and examples from repository.
- 2022-07-28: Add Docker recipe used for CI.
- 2022-07-27: Add DuMuX solverdummy and add it as test case.
- 2022-07-27: Add support for exchanging vector quantities.
- 2022-07-27: Make sure clang-format fails when files a badly formatted.
- 2022-05-25: Add CMake guards to prevent build targe generation of cases that depend of `dune-subgrid`, if `dune-subgrid` is not installed.
- 2022-05-24: Added missing include of `limits` in `couplingadapter.cc`.
- 2022-05-17: Added base setup for extended documentation to be hosted on ReadTheDocs and being created by `mkdocs`. Also adds a base configuration and style of the documentation.
- 2022-05-17: Added configuation for Markdown linter `markdownlint` and added it to CI. The linter can be called locally by typing `mdl .` from the root of the repository. This also led to an updated configuration of the CI.
- 2022-05-17: Add some more documentation on how to install the adapter to the `README.md`.
- 2022-02-18: Updated CI to use images from account `ajaust` from Dockerhub. Changed tolerance for partitioned tests to 5e-5 due to minimal changes in the solution with the new images on a new VM.
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
