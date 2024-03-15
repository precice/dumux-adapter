add_executable(ff_flow_over_cube_3d EXCLUDE_FROM_ALL main_ff-reversed.cc)
add_executable(pm_flow_over_cube_3d EXCLUDE_FROM_ALL main_pm-reversed.cc)

# add a symlink for each input file
add_input_file_links()
# add a symlink for each preCICE configuration file
add_precice_file_links()

# Dummy test is needed to build porous media solver. The default
# `dumux_add_test` (and the underlying `dune_add_test`) only accept
# one target executable.
dumux_add_test(NAME dummytarget_to_build_pm_flow_over_cube_3d
              TARGET pm_flow_over_cube_3d
              LABELS freeflow stokes precice darcy 3d
              TIMEOUT 5
              COMMAND ${CMAKE_SOURCE_DIR}/test/return-test-passed.sh
)


dumux_add_test(NAME test_ff_pm_flow_over_cube_3d
              TARGET ff_flow_over_cube_3d
              LABELS freeflow stokes precice darcy
              TIMEOUT 30
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-test.py
              CMD_ARGS
              --solver-one ${CMAKE_CURRENT_BINARY_DIR}/ff_flow_over_cube_3d
              --solver-two ${CMAKE_CURRENT_BINARY_DIR}/pm_flow_over_cube_3d
              --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-cube-3d/test_stokes_first_si_stokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-3d-precice-stokes-00005.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-cube-3d/test_stokes_first_si_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-3d-precice-darcy-00005.vtu
              --dumux-param-file ${CMAKE_CURRENT_BINARY_DIR}/params.input
              --precice-config-file ${CMAKE_CURRENT_BINARY_DIR}/precice-config.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-cube-3d/test_stokes_first_si_coupliter_stokes.log ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-cube_3d"
              --relative 5e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_flow_over_cube_3d PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/test:$ENV{PYTHONPATH})
