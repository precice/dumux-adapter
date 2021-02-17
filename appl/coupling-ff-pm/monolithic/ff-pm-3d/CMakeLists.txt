
add_input_file_links()

add_executable(ff-pm-3d EXCLUDE_FROM_ALL fvca-monolithic-3d.cc)

dune_add_test(NAME test_ff-pm-3d-stokes
              TARGET ff-pm-3d
              LABELS ff-pm stokes monolithic darcy
              TIMEOUT 10
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/runtest.py
              CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/monolithic/ff-pm-3d/ref-ff-pm-3d-stokes_stokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-3d-stokes_stokes-00001.vtu
                      ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/monolithic/ff-pm-3d/ref-ff-pm-3d-stokes_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-3d-stokes_darcy-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/ff-pm-3d ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/monolithic/ff-pm-3d/ff-pm-3d-stokes.input"
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-8,"p":1e-8}
)
set_tests_properties(test_ff-pm-3d-stokes PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})


dune_add_test(NAME test_ff-pm-3d-navier-stokes
              TARGET ff-pm-3d
              LABELS ff-pm stokes monolithic darcy
              TIMEOUT 10
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/runtest.py
              CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/monolithic/ff-pm-3d/ref-ff-pm-3d-navier-stokes_stokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-3d-navier-stokes_stokes-00001.vtu
                      ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/monolithic/ff-pm-3d/ref-ff-pm-3d-navier-stokes_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-3d-navier-stokes_darcy-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/ff-pm-3d ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/monolithic/ff-pm-3d/ff-pm-3d-navier-stokes.input"
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-8,"p":1e-8}
)
set_tests_properties(test_ff-pm-3d-navier-stokes PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})