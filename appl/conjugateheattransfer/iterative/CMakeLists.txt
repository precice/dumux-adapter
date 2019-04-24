add_input_file_links()
dune_symlink_to_source_files(FILES "precice-config.xml" "precice-config-serial-implicit.xml")

add_executable(test_freeflow EXCLUDE_FROM_ALL main_freeflow.cc precicewrapper.cc dumuxpreciceindexwrapper.cc)
add_executable(test_solidenergy EXCLUDE_FROM_ALL main_solidenergy.cc precicewrapper.cc dumuxpreciceindexwrapper.cc)
#add_executable(test_precicewrapper EXCLUDE_FROM_ALL precicewrapper.cc)

target_compile_definitions(test_freeflow PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(test_solidenergy PUBLIC "ENABLEMONOLITHIC=0")

#set_property( TARGET libprecice IMPORTED_LOCATION $ENV{PRECICE_ROOT}/build/last/librprecice.so )

find_library(PRECICE_LIB precice)
target_link_libraries(test_freeflow precice)
target_link_libraries(test_solidenergy precice)
#target_link_libraries( test_freeflow precice )
#set_property( TARGET libprecice IMPORTED_LOCATION $ENV{PRECICE_ROOT}/build/last/librprecice.so )
