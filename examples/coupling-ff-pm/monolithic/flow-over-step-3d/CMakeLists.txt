add_executable(flow-over-step-3d EXCLUDE_FROM_ALL main.cc)

# add a symlink for each input file
add_input_file_links()


#dune_add_test(NAME test_flow-over-step-3d
#              TARGET flow-over-step-3d
#              LABELS ff-pm stokes monolithic darcy
#              TIMEOUT 10
#              CMAKE_GUARD HAVE_UMFPACK
#              COMMAND ${CMAKE_SOURCE_DIR}/test/runtest.py
#              CMD_ARGS --script fuzzy
#                      --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/monolithic/flow-over-step-3d/ref-sol-stokes-flow-over-step-3d.vtu
#                      ${CMAKE_CURRENT_BINARY_DIR}/flow-over-step-3d_stokes-00001.vtu
#                      ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/monolithic/flow-over-step-3d/ref-sol-darcy-flow-over-step-3d.vtu
#                      ${CMAKE_CURRENT_BINARY_DIR}/flow-over-step-3d_darcy-00001.vtu
#              --command "${CMAKE_CURRENT_BINARY_DIR}/flow-over-step-3d ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/monolithic/flow-over-step-3d/ref-flow-over-step-3d.input"
#              --relative 1e-6
#              --zeroThreshold {"velocity_liq \(m/s\)":1e-6,"p":1e-8}
#)
#set_tests_properties(test_flow-over-step-3d PROPERTIES
#    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})
