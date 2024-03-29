project(dumux-precice CXX)
# We require version CMake version 3.1 to prevent issues
# with dune_enable_all_packages and older CMake versions.
cmake_minimum_required(VERSION 3.13)

if(NOT (dune-common_DIR
        OR dune-common_ROOT
        OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME}
      dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
# Find preCICE library
find_package(precice 3.0.0 REQUIRED CONFIG)

list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")
#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

dune_enable_all_packages(MODULE_LIBRARIES dumux-precice)

# enforce C++-17
dune_require_cxx_standard(MODULE "DuMuX-preCICE" VERSION 17)

add_subdirectory(cmake/modules)
add_subdirectory(doc)
add_subdirectory(dumux-precice)
add_subdirectory(examples)
add_subdirectory(test EXCLUDE_FROM_ALL)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
