#add_input_file_links()
dune_symlink_to_source_files(FILES "../monolithic/params.input")

add_executable(test_conjugate_heat_monolithic_const_dt EXCLUDE_FROM_ALL main.cc)
