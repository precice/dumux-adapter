dune_add_test(NAME orig_turbulence_coupling_ff-pm
              SOURCES ex_turbulence_coupling_ff-pm.cc
              COMPILE_DEFINITIONS EXNUMBER=0)

dune_add_test(NAME sol_a_ex_turbulence_coupling_ff-pm
              SOURCES ex_turbulence_coupling_ff-pm.cc
              COMPILE_DEFINITIONS EXNUMBER=1)

dune_add_test(NAME sol_b_ex_turbulence_coupling_ff-pm
              SOURCES ex_turbulence_coupling_ff-pm.cc
              COMPILE_DEFINITIONS EXNUMBER=2)

# add exercise to the common target
add_dependencies(test_exercises sol_a_ex_turbulence_coupling_ff-pm sol_b_ex_turbulence_coupling_ff-pm sol_b_ex_turbulence_coupling_ff-pm)

# add a symlink for each input file
add_input_file_links()
