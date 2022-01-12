
add_input_file_links()

add_executable(test_ff_pm_flow_over_square_2d EXCLUDE_FROM_ALL main.cc)

dune_add_test(NAME test_ff_pm_mono_flow_over_square_2d
              TARGET test_ff_pm_flow_over_square_2d
              LABELS ff-pm stokes monolithic darcy
              TIMEOUT 10
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/runtest.py
              CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/monolithic/flow-over-square-2d/stokes-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-darcy-monolithic-coupling-test_stokes-00001.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/monolithic/flow-over-square-2d/darcy-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-darcy-monolithic-coupling-test_darcy-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_pm_flow_over_square_2d ${CMAKE_CURRENT_BINARY_DIR}/params.input"
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_mono_flow_over_square_2d PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})