# executables for the monolithic case
dune_add_test(NAME ff-pm_monolithic
              SOURCES main.cc)

# add a symlink for each input file
add_input_file_links()
