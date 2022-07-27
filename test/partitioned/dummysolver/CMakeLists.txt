add_executable(test_partitioned_dummysolver EXCLUDE_FROM_ALL main_dummysolver.cc)

target_link_libraries(test_partitioned_dummysolver PRIVATE dumux-precice)

dune_symlink_to_source_files(FILES precice-dummy-solver-config.xml)
dune_symlink_to_source_files(FILES Allrun.sh)

dumux_add_test(NAME test_partitioned_dummysolver
              TARGET test_partitioned_dummysolver
              LABELS dummy precice partitioned
              TIMEOUT 30
              COMMAND test_partitioned_dummysolver)
