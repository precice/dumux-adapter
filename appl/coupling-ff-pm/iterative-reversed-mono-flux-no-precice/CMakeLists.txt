# executables for the iterative case

add_executable(test_ff_reversed_mono_flux-no-precice  EXCLUDE_FROM_ALL main_ff-mono-flux.cc )
add_executable(test_pm_reversed_mono_flux-no-precice  EXCLUDE_FROM_ALL main_pm-mono-flux.cc )

target_compile_definitions(test_ff_reversed_mono_flux-no-precice PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(test_pm_reversed_mono_flux-no-precice  PUBLIC "ENABLEMONOLITHIC=0")

# add a symlink for each input file
add_input_file_links()

macro(add_precice_file_links)
  FILE(GLOB precice_files *.xml)
  foreach(VAR ${input_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
add_precice_file_links()
