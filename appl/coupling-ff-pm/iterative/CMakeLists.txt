# executables for the iterative case
dune_add_test(NAME test_ff
              SOURCES main_ff.cc
              COMPILE_DEFINITIONS ENABLEMONOLITHIC=0)

dune_add_test(NAME test_pm
              SOURCES main_pm.cc
              COMPILE_DEFINITIONS ENABLEMONOLITHIC=0)

# add a symlink for each input file
add_input_file_links()
