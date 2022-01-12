#dune_add_test(NAME exercise_interface_coupling_ff-pm_original
#              SOURCES main.cc
#              COMPILE_DEFINITIONS EXNUMBER=0)

#dune_add_test(NAME exercise_interface_coupling_ff-pm_a_solution
#              SOURCES main.cc
#              COMPILE_DEFINITIONS EXNUMBER=1)

#dune_add_test(NAME exercise_interface_coupling_ff-pm_b_solution
#              SOURCES main.cc
#              COMPILE_DEFINITIONS EXNUMBER=2)

#dune_add_test(NAME flow-over-step-2d
#              SOURCES main.cc)


# add exercise to the common target
#add_dependencies(test_exercises exercise_interface_coupling_ff-pm_a_solution exercise_interface_coupling_ff-pm_b_solution exercise_interface_coupling_ff-pm_c_solution)
#add_dependencies(flow-over-step-2d)
add_executable(flow-over-step-2d EXCLUDE_FROM_ALL main.cc)

# add a symlink for each input file
add_input_file_links()


dune_add_test(NAME test_flow-over-step-2d
              TARGET flow-over-step-2d
              LABELS ff-pm stokes monolithic darcy
              TIMEOUT 10
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/runtest.py
              CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/monolithic/flow-over-step-2d/ref-sol-stokes-flow-over-step-2d.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/flow-over-step-2d_stokes-00001.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/monolithic/flow-over-step-2d/ref-sol-darcy-flow-over-step-2d.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/flow-over-step-2d_darcy-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/flow-over-step-2d ${CMAKE_SOURCE_DIR}/test/reference-solutions/monolithic/flow-over-step-2d/ref-flow-over-step-2d.input"
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-6,"p":1e-8}
)
set_tests_properties(test_flow-over-step-2d PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})
