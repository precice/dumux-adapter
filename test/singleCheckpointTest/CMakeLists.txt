add_executable(dumuxprecice_singlecheckpointtest EXCLUDE_FROM_ALL main.cc)

dune_symlink_to_source_files(FILES test.xml Allrun.sh params_one.input params_two.input)

dumux_add_test(NAME dumuxprecice_singlecheckpointtest
              TARGET dumuxprecice_singlecheckpointtest
              LABELS singlecheckpointtest precice dummy
              TIMEOUT 30
              COMMAND Allrun.sh)
