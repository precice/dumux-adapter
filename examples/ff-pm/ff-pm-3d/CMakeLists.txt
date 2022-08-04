add_executable(ff_flow_over_cube_3d EXCLUDE_FROM_ALL main_ff-reversed.cc)
add_executable(pm_flow_over_cube_3d EXCLUDE_FROM_ALL main_pm-reversed.cc)

target_link_libraries(ff_flow_over_cube_3d PRIVATE dumux-precice)
target_link_libraries(pm_flow_over_cube_3d PRIVATE dumux-precice)

# add a symlink for each input file
add_input_file_links()
# add a symlink for each preCICE configuration file
add_precice_file_links()

# Dummy test is needed to build porous media solver. The default
# `dumux_add_test` (and the underlying `dune_add_test`) only accept
# one target executable.
dumux_add_test(NAME dummytarget_to_build_pm_flow_over_cube_3d
              TARGET pm_flow_over_cube_3d
              LABELS freeflow stokes precice darcy partitioned 3d
              TIMEOUT 5
              COMMAND ${CMAKE_SOURCE_DIR}/test/return-test-passed.sh
)

#dumux_add_test(NAME dummytarget_to_build_ff_flow_over_cube_3d
#              TARGET ff_flow_over_cube_3d
#              LABELS freeflow stokes precice darcy partitioned 3d
#              TIMEOUT 5
#              COMMAND ${CMAKE_SOURCE_DIR}/test/return-test-passed.sh
#)


dumux_add_test(NAME test_ff_pm_flow_over_cube_3d
              TARGET ff_flow_over_cube_3d
              TARGET pm_flow_over_cube_3d
              LABELS freeflow stokes precice darcy partitioned
              TIMEOUT 30
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-test.py
              CMD_ARGS
              --skip-checks
              --solver-one ${CMAKE_CURRENT_BINARY_DIR}/ff_flow_over_cube_3d
              --solver-two ${CMAKE_CURRENT_BINARY_DIR}/pm_flow_over_cube_3d
              --files DUMMY_VTU1
                      DUMMY_VTU2
                      DUMMY_VTU3
                      DUMMY_VTU4
              --dumux-param-file ${CMAKE_CURRENT_BINARY_DIR}/params-stokesdarcy.input
              --precice-config-file ${CMAKE_CURRENT_BINARY_DIR}/precice-config.xml
              --precice-iteration-files DUMMY_PRECICE_ITERATION1 DUMMY_PRECICE_ITERATION2
              --case-name "flow-over-cube_3d"
              --relative 5e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_flow_over_cube_3d PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/test:$ENV{PYTHONPATH})