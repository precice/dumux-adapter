add_executable(test_ff_reversed-3d EXCLUDE_FROM_ALL main_ff-reversed.cc)
add_executable(test_pm_reversed-3d EXCLUDE_FROM_ALL main_pm-reversed.cc)

target_compile_definitions(test_ff_reversed-3d PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(test_pm_reversed-3d PUBLIC "ENABLEMONOLITHIC=0")

target_link_libraries(test_ff_reversed-3d PRIVATE dumux-precice)
target_link_libraries(test_pm_reversed-3d PRIVATE dumux-precice)

# add a symlink for each input file
add_input_file_links()
# add a symlink for each preCICE configuration file
add_precice_file_links()

dune_symlink_to_source_files(FILES run-iterative-simulation-test.sh)
dune_symlink_to_source_files(FILES precice-config.xml)
