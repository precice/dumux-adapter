add_executable(ff_flow_over_cube_3d EXCLUDE_FROM_ALL main_ff-reversed.cc)
add_executable(pm_flow_over_cube_3d EXCLUDE_FROM_ALL main_pm-reversed.cc)

target_link_libraries(ff_flow_over_cube_3d PRIVATE dumux-precice)
target_link_libraries(pm_flow_over_cube_3d PRIVATE dumux-precice)

# add a symlink for each input file
add_input_file_links()
# add a symlink for each preCICE configuration file
add_precice_file_links()

dune_symlink_to_source_files(FILES run-iterative-simulation-test.sh)
dune_symlink_to_source_files(FILES precice-config.xml)

# Dummy test is needed to build porous media solver. The default
# `dumux_add_test` (and the underlying `dune_add_test`) only accept
# one target executable.
dumux_add_test(NAME dummy_test_ff_pm_part_flow_over_square_2d
              TARGET pm_flow_over_square_2d
              LABELS freeflow stokes precice darcy partitioned
              TIMEOUT 30
              COMMAND ${CMAKE_SOURCE_DIR}/test/return-test-passed.sh
)