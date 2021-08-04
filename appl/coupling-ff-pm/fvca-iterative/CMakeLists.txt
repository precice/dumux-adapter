# executables for the iterative case
set(LIBRARY_DIRS "")
if(DEFINED ENV{LD_LIBRARY_PATH})
    message(STATUS "LD_LIBRARY_PATH is set. Use it for hinting")
    string(REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH} )
endif()
#message(STATUS "LIBRARY_DIRS: ${LIBRARY_DIRS}")
find_library(PRECICE_LIB NAMES precice HINTS ${LIBRARY_DIRS})
#find_package(Boost 1.65.1 REQUIRED COMPONENTS log system) #Require same version as preCICE


add_executable(fvca-iterative-ff EXCLUDE_FROM_ALL main_fvca-iterative-ff.cc ../../precice-adapter/src/preciceadapter.cc ../../precice-adapter/src/dumuxpreciceindexwrapper.cc)
add_executable(fvca-iterative-pm EXCLUDE_FROM_ALL main_fvca-iterative-pm.cc ../../precice-adapter/src/preciceadapter.cc ../../precice-adapter/src/dumuxpreciceindexwrapper.cc)

target_compile_definitions(fvca-iterative-ff PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(fvca-iterative-pm PUBLIC "ENABLEMONOLITHIC=0")

target_link_libraries(fvca-iterative-ff PRIVATE ${PRECICE_LIB})
target_link_libraries(fvca-iterative-pm PRIVATE ${PRECICE_LIB})
add_input_file_links()

macro(add_precice_file_links)
  FILE(GLOB precice_files *.xml)
  foreach(VAR ${input_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
add_precice_file_links()


dune_symlink_to_source_files(FILES precice-config-serial-implicit-stokes-first.xml)
dune_symlink_to_source_files(FILES precice-config-serial-implicit-darcy-first.xml)
dune_symlink_to_source_files(FILES precice-config-parallel-implicit.xml)
dune_symlink_to_source_files(FILES Allrun.sh)
dune_symlink_to_source_files(FILES Allclean.sh)


dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_stokes_first
              TARGET fvca-iterative-ff
              TARGET fvca-iterative-pm
              LABELS freeflow stokes precice darcy
              TIMEOUT 15
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/serial-implicit-stokes-first/stokes-iterative-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-iterative-00005.vtu
                      ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/serial-implicit-stokes-first/darcy-iterative-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00005.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/params.input
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/serial-implicit-stokes-first/precice-config.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/serial-implicit-stokes-first/precice-FreeFlow-iterations.log ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-box-2d-si-stokes-first"
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_serial_implicit_stokes_first PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_darcy_first
    TARGET fvca-iterative-ff
    TARGET fvca-iterative-pm
    LABELS freeflow darcy precice darcy
    TIMEOUT 15
    CMAKE_GUARD HAVE_UMFPACK
    COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py
    CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/serial-implicit-darcy-first/stokes-iterative-final.vtu
            ${CMAKE_CURRENT_BINARY_DIR}/stokes-iterative-00006.vtu
            ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/serial-implicit-darcy-first/darcy-iterative-final.vtu
            ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00006.vtu
    --dumux-param-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/params.input
    --precice-config-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/serial-implicit-darcy-first/precice-config.xml
    --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/serial-implicit-darcy-first/precice-Darcy-iterations.log
        ${CMAKE_CURRENT_BINARY_DIR}/precice-Darcy-iterations.log
    --case-name "flow-over-box-2d-si-darcy-first"
    --relative 1e-6
    --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_serial_implicit_darcy_first PROPERTIES
ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_coupling_parallel_implicit
              TARGET fvca-iterative-ff
              TARGET fvca-iterative-pm
              LABELS freeflow stokes precice darcy
              TIMEOUT 15
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/parallel-implicit/stokes-iterative-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-iterative-00007.vtu
                      ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/parallel-implicit/darcy-iterative-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00007.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/params.input
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/parallel-implicit/precice-config.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/partitioned/2d/flow-over-box/parallel-implicit/precice-FreeFlow-iterations.log
              ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-box-2d-pi"
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
# Adding DuMuX directory to Python search path such that we can use DuMuX'x fuzzy (inexact) testing of vtu files
set_tests_properties(test_ff_pm_coupling_parallel_implicit PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

    ### Setting test dependencies
set_property(TEST test_ff_pm_coupling_serial_implicit_darcy_first APPEND PROPERTY DEPENDS test_ff_pm_coupling_serial_implicit_stokes_first)
set_property(TEST test_ff_pm_coupling_parallel_implicit APPEND PROPERTY DEPENDS test_ff_pm_coupling_serial_implicit_darcy_first)
