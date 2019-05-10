# executables for ex_interface_coupling_ff-pm
dune_add_test(NAME ex_turbulence_coupling_ff-pm
              SOURCES ex_turbulence_coupling_ff-pm.cc)


# add a symlink for each input file
add_input_file_links()

# add tutorial to the common target
add_dependencies(test_exercises ex_turbulence_coupling_ff-pm)
