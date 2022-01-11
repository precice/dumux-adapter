add_executable(fvca-iterative-ff EXCLUDE_FROM_ALL main_fvca-iterative-ff.cc)
add_executable(fvca-iterative-pm EXCLUDE_FROM_ALL main_fvca-iterative-pm.cc)

target_compile_definitions(fvca-iterative-ff PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(fvca-iterative-pm PUBLIC "ENABLEMONOLITHIC=0")

target_link_libraries(fvca-iterative-ff PRIVATE dumux-precice)
target_link_libraries(fvca-iterative-pm PRIVATE dumux-precice)

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
              TIMEOUT 30
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_stokes_first_si_stokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-iterative-00005.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_stokes_first_si_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00005.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/params-stokesdarcy.input
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/precice-config-si-free-flow-first.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_stokes_first_si_coupliter_stokes.log ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-box-2d-si-stokes-first"
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_serial_implicit_stokes_first PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_stokes_second
    TARGET fvca-iterative-ff
    TARGET fvca-iterative-pm
    LABELS freeflow darcy precice darcy
    TIMEOUT 30
    CMAKE_GUARD HAVE_UMFPACK
    COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py
    CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_stokes_second_si_stokes.vtu
              ${CMAKE_CURRENT_BINARY_DIR}/stokes-iterative-00006.vtu
              ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_stokes_second_si_darcy.vtu
              ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00006.vtu
            --dumux-param-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/params-stokesdarcy.input
            --precice-config-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/precice-config-si-free-flow-second.xml
            --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_stokes_second_si_coupliter_darcy.log ${CMAKE_CURRENT_BINARY_DIR}/precice-Darcy-iterations.log
    --case-name "flow-over-box-2d-si-stokes-second"
    --relative 1e-6
    --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_serial_implicit_stokes_second PROPERTIES
ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_coupling_parallel_implicit_stokes
              TARGET fvca-iterative-ff
              TARGET fvca-iterative-pm
              LABELS freeflow stokes precice darcy
              TIMEOUT 30
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_stokes_pi_stokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-iterative-00007.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_stokes_pi_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00007.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/params-stokesdarcy.input
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/precice-config-pi.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_stokes_pi_coupliter_stokes.log
              ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-box-2d-pi"
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
# Adding DuMuX directory to Python search path such that we can use DuMuX'x fuzzy (inexact) testing of vtu files
set_tests_properties(test_ff_pm_coupling_parallel_implicit_stokes PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

### Setting test dependencies
set_property(TEST test_ff_pm_coupling_serial_implicit_stokes_second APPEND PROPERTY DEPENDS test_ff_pm_coupling_serial_implicit_stokes_first)
set_property(TEST test_ff_pm_coupling_parallel_implicit_stokes APPEND PROPERTY DEPENDS test_ff_pm_coupling_serial_implicit_stokes_second)

######################
# Navier-Stokes tests
######################

dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_navierstokes_first
              TARGET fvca-iterative-ff
              TARGET fvca-iterative-pm
              LABELS freeflow navierstokes precice darcy
              TIMEOUT 60
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_navierstokes_first_si_navierstokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/navierstokes-iterative-00006.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_navierstokes_first_si_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00006.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/params-navierstokesdarcy.input
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/precice-config-si-free-flow-first.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_navierstokes_first_si_coupliter_navierstokes.log ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-box-2d-si-navierstokes-first"
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_serial_implicit_navierstokes_first PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_navierstokes_second
    TARGET fvca-iterative-ff
    TARGET fvca-iterative-pm
    LABELS freeflow darcy precice darcy
    TIMEOUT 60
    CMAKE_GUARD HAVE_UMFPACK
    COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py
    CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_navierstokes_second_si_navierstokes.vtu
              ${CMAKE_CURRENT_BINARY_DIR}/navierstokes-iterative-00007.vtu
              ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_navierstokes_second_si_darcy.vtu
              ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00007.vtu
            --dumux-param-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/params-navierstokesdarcy.input
            --precice-config-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/precice-config-si-free-flow-second.xml
            --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_navierstokes_second_si_coupliter_darcy.log ${CMAKE_CURRENT_BINARY_DIR}/precice-Darcy-iterations.log
    --case-name "flow-over-box-2d-si-navierstokes-second"
    --relative 1e-6
    --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_serial_implicit_navierstokes_second PROPERTIES
ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_coupling_parallel_implicit_navierstokes
              TARGET fvca-iterative-ff
              TARGET fvca-iterative-pm
              LABELS freeflow navierstokes precice darcy
              TIMEOUT 60
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_navierstokes_pi_navierstokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/navierstokes-iterative-00011.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_navierstokes_pi_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00011.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/params-navierstokesdarcy.input
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/partitioned/flow-over-square-2d/precice-config-pi.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/partitioned/flow-over-square-2d/test_navierstokes_pi_coupliter_navierstokes.log
              ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-box-2d-pi"
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
# Adding DuMuX directory to Python search path such that we can use DuMuX'x fuzzy (inexact) testing of vtu files
set_tests_properties(test_ff_pm_coupling_parallel_implicit_navierstokes PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

# Adding DuMuX directory to Python search path such that we can use DuMuX'x fuzzy (inexact) testing of vtu files
set_property(TEST test_ff_pm_coupling_serial_implicit_navierstokes_first APPEND PROPERTY DEPENDS test_ff_pm_coupling_parallel_implicit_stokes)
set_property(TEST test_ff_pm_coupling_serial_implicit_navierstokes_second APPEND PROPERTY DEPENDS test_ff_pm_coupling_serial_implicit_navierstokes_first)
set_property(TEST test_ff_pm_coupling_parallel_implicit_navierstokes APPEND PROPERTY DEPENDS test_ff_pm_coupling_serial_implicit_navierstokes_second)