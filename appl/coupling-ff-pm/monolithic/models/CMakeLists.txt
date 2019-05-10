# executables for ex_interface_coupling_ff-pm
dune_add_test(NAME ex_models_coupling_ff-pm
              SOURCES ex_models_coupling_ff-pm.cc)

# add tutorial to the common target
add_dependencies(test_exercises ex_models_coupling_ff-pm)

# add a symlink for each input file
add_input_file_links()
