dune_add_test(NAME orig_ex_interface_coupling_ff-pm
              SOURCES ex_interface_coupling_ff-pm.cc
              COMPILE_DEFINITIONS EXNUMBER=0)

dune_add_test(NAME sol_a_ex_interface_coupling_ff-pm
              SOURCES ex_interface_coupling_ff-pm.cc
              COMPILE_DEFINITIONS EXNUMBER=1)

dune_add_test(NAME sol_b_ex_interface_coupling_ff-pm
              SOURCES ex_interface_coupling_ff-pm.cc
              COMPILE_DEFINITIONS EXNUMBER=2)

dune_add_test(NAME sol_c_ex_interface_coupling_ff-pm
              SOURCES ex_interface_coupling_ff-pm.cc
              COMPILE_DEFINITIONS EXNUMBER=3)


# add exercise to the common target
add_dependencies(test_exercises sol_a_ex_interface_coupling_ff-pm sol_b_ex_interface_coupling_ff-pm sol_b_ex_interface_coupling_ff-pm)

# add a symlink for each input file
add_input_file_links()
