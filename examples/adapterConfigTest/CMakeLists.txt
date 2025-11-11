add_executable(adapterConfigTest EXCLUDE_FROM_ALL main.cc)

dune_symlink_to_source_files(FILES precice-config.xml)
dune_symlink_to_source_files(FILES Allrun.sh)

dumux_add_test(NAME adapterConfigTest
              TARGET adapterConfigTest
              LABELS test precice
              TIMEOUT 30
              COMMAND Allrun.sh)
