add_input_file_links()

add_executable(test_freeflow EXCLUDE_FROM_ALL main_freeflow.cc)
add_executable(test_solidenergy EXCLUDE_FROM_ALL main_solidenergy.cc)

target_compile_definitions(test_freeflow PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(test_solidenergy PUBLIC "ENABLEMONOLITHIC=0")

#set_property( TARGET libprecice IMPORTED_LOCATION $ENV{PRECICE_ROOT}/build/last/librprecice.so )

find_library(PRECICE_LIB precice)
target_link_libraries(test_freeflow precice)
target_link_libraries(test_solidenergy precice)
#target_link_libraries( test_freeflow precice )
#set_property( TARGET libprecice IMPORTED_LOCATION $ENV{PRECICE_ROOT}/build/last/librprecice.so )
