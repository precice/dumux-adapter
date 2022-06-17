# Create file links to preCICE configuration files
macro(add_precice_file_links)
  FILE(GLOB precice_input_files *.xml)
  foreach(VAR ${precice_input_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
