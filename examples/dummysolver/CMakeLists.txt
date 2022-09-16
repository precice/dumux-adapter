add_executable(dumuxprecice_dummysolver EXCLUDE_FROM_ALL main_dummysolver.cc)

target_link_libraries(dumuxprecice_dummysolver PRIVATE dumux-precice)

dune_symlink_to_source_files(FILES precice-dummy-solver-config.xml)
dune_symlink_to_source_files(FILES Allrun.sh)

dumux_add_test(NAME dumuxprecice_dummysolver
              TARGET dumuxprecice_dummysolver
              LABELS dummy precice
              TIMEOUT 30
              COMMAND Allrun.sh)
